# -*- coding: utf-8 -*-
"""The process pipeline the FASTQ and FAST5 converters share.

A lister process queues work items, a pool of parser processes turns each
item into one record per read, and a writer process collects the records
into the summary file:

    lister --in_q--> parser x N --out_q--> writer
                \\________ error_q ________/  -> parent

Each converter supplies only what differs: the work items, how an item
becomes records, and any last change to the table before it is written.
Starting the processes, passing errors back, counting fields and stopping
everything on the way out live here, once.
"""

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import logging
from collections import Counter
from time import time

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from pycoQC.common import (
    dict_to_str, format_worker_error, get_logger, pycoQCError, stop_processes,
    wait_for_workers, worker_context)


def _new_counts():
    return {"overall": Counter(), "fields_found": Counter(),
            "fields_not_found": Counter()}


#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class ConverterPipeline:
    """Base class for converters that parse their input in worker processes.

    A subclass sets, before calling `_run`:
        seq_summary_fn, verbose_level, log_level
        threads: the number of parser processes
    and implements:
        _work_items(): yields picklable work items, in the lister
        _records(item, counts): yields one dict per read, in a parser,
            counting fields in counts["fields_found"] and
            counts["fields_not_found"]
        _finalise(df): optional, the table just before it is written

    The stage methods run in other processes, so they reach the converter
    as a pickled copy: anything they need must be an instance attribute set
    before `_run`, not state changed on the class in the parent.
    """

    # What the input is called in messages: "fastq", "fast5".
    kind = None

    # Work items waiting for a parser. A large item (a block of lines) needs
    # a short queue to bound memory; a file path does not.
    in_queue_size = 1000

    # Records travel to the writer in lists of up to this many. Sending each
    # read alone spent more time pickling than the writer spent on the rest.
    RECORDS_PER_BATCH = 1000

    #~~~~~~~~~~~~~~TO IMPLEMENT~~~~~~~~~~~~~~#
    def _work_items(self):
        raise NotImplementedError

    def _records(self, item, counts):
        raise NotImplementedError

    def _finalise(self, df):
        return df

    #~~~~~~~~~~~~~~PARENT~~~~~~~~~~~~~~#
    @property
    def _log(self):
        # Each converter logs under its own module's name, as it always has.
        return logging.getLogger(type(self).__module__)

    def _run(self):
        """Start the lister, parsers and writer, and wait for them to finish."""
        ctx = worker_context()
        in_q = ctx.Queue(maxsize=self.in_queue_size)
        out_q = ctx.Queue(maxsize=max(8, 4 * self.threads))
        error_q = ctx.Queue()
        counter_q = ctx.Queue()

        ps_list = [ctx.Process(target=self._lister, args=(in_q, error_q))]
        for i in range(self.threads):
            ps_list.append(ctx.Process(
                target=self._parser, args=(in_q, out_q, error_q, counter_q, i)))
        ps_list.append(ctx.Process(
            target=self._writer, args=(out_q, error_q, counter_q)))

        self._log.info("Start processing {} files".format(self.kind))
        try:
            for proc in ps_list:
                proc.start()
            wait_for_workers(ps_list, error_q)
            for proc in ps_list:
                proc.join()

        except BaseException:
            self._log.info("\nAn error occurred. All processes were killed\n")
            raise

        finally:
            # Every exit from here, including a clean one, leaves nothing
            # running. Listing the exceptions to clean up after meant any
            # other kind left the whole pool alive.
            stop_processes(ps_list)

    #~~~~~~~~~~~~~~WORKERS~~~~~~~~~~~~~~#
    def _start_worker(self):
        # A worker starts with fresh loggers rather than the parent's.
        return get_logger(type(self).__module__, level=self.log_level)

    def _lister(self, in_q, error_q):
        """Queue the work items, then one stop signal per parser."""
        log = self._start_worker()
        try:
            n = 0
            for item in self._work_items():
                in_q.put(item)
                n += 1
            log.debug("[LISTER] Queued {} work items".format(n))
        except Exception as E:
            error_q.put(format_worker_error(E))
        finally:
            for _ in range(self.threads):
                in_q.put(None)

    def _parser(self, in_q, out_q, error_q, counter_q, worker_id):
        """Turn work items into batches of records for the writer."""
        log = self._start_worker()
        log.debug("[PARSER_{:02}] Start".format(worker_id))
        try:
            counts = _new_counts()
            for item in iter(in_q.get, None):
                batch = []
                for record in self._records(item, counts):
                    if not record:
                        counts["overall"]["invalid reads"] += 1
                        continue
                    counts["overall"]["valid reads"] += 1
                    batch.append(record)
                    if len(batch) >= self.RECORDS_PER_BATCH:
                        out_q.put(batch)
                        batch = []
                if batch:
                    out_q.put(batch)
            counter_q.put(counts)
        except Exception as E:
            error_q.put(format_worker_error(E))
        finally:
            out_q.put(None)

    def _writer(self, out_q, error_q, counter_q):
        """Collect every record, write the summary and report the counts.

        Signals the parent with None on error_q once done, after the error
        if there was one.
        """
        log = self._start_worker()
        t = time()
        try:
            rows = []
            with tqdm(unit=" reads", mininterval=0.1, smoothing=0.1,
                      disable=self.verbose_level == 2) as pbar:
                for _ in range(self.threads):
                    for batch in iter(out_q.get, None):
                        rows.extend(batch)
                        pbar.update(len(batch))

            # Files holding no reads used to give an empty summary, which
            # pycoQC then rejected with a message about the summary rather
            # than about the input files.
            log.debug("[WRITER] Write data to file")
            if not rows:
                raise pycoQCError(
                    "No reads found in the {0} file(s). The files were read "
                    "but hold no reads.".format(self.kind))
            df = self._finalise(pd.DataFrame(rows))
            df.to_csv(self.seq_summary_fn, sep="\t", index=False)

            counts = _new_counts()
            for _ in range(self.threads):
                for key, counter in counter_q.get().items():
                    counts[key].update(counter)
            log.info("Overall counts {}".format(dict_to_str(counts["overall"])))
            log.info("fields found {}".format(dict_to_str(counts["fields_found"])))
            log.info("fields not found {}".format(dict_to_str(counts["fields_not_found"])))

            log.warning("Total reads: {} / Average speed: {} reads/s\n".format(
                len(df), round(len(df) / (time() - t), 2)))

        except Exception as E:
            error_q.put(format_worker_error(E))
        finally:
            error_q.put(None)

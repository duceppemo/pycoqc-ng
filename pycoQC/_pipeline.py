# -*- coding: utf-8 -*-
"""The process pipeline every converter shares.

A lister process queues work items, a pool of parser processes turns each
item into one record per read, and a writer process collects the records
into the summary file:

    lister --in_q--> parser x N --out_q--> writer
                \\________ error_q ________/  -> parent

Each converter supplies only what differs: the work items, how an item
becomes records, and any last change to the table before it is written.
Starting the processes, passing errors back, putting the rows in a fixed
order, counting fields, reporting and stopping everything on the way out
live here, once.

Everything the user is told is logged by the parent, from a report the
writer sends back, so it reaches the caller's logging configuration rather
than only a worker's stderr.
"""

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import logging
import multiprocessing as mp
import os
import queue
import shutil
import sys
import tempfile
import threading
from collections import Counter
from contextlib import contextmanager
from time import time

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from pycoQC._worker_env import WORKER_ENV
from pycoQC.common import (
    VERBOSE_LEVELS, check_converter_io, dict_to_str, format_worker_error,
    get_logger, pycoQCError, recursive_file_gen)


#~~~~~~~~~~~~~~PROCESS HELPERS~~~~~~~~~~~~~~#
def worker_context():
    """The multiprocessing context the converters start their workers in.

    Forking copies the parent's threads' locks in whatever state they are
    in, and a child can deadlock on one another thread held at the time.
    Python 3.12 warns about it whenever the parent has a thread running, as
    it does under pytest, in a notebook or once numpy or pandas has started
    one. A fork server forks from a clean single threaded process instead,
    and is the Linux default from Python 3.14.

    Elsewhere, spawn, as Python itself does. Windows has no fork server, and
    on macOS forking a process that has loaded the system frameworks is
    unsafe whatever its threads. The fork server still held a thread there,
    BLAS pools capped or not, and warned on every fork.

    The server preloads pycoQC through _forkserver_preload, which keeps
    numpy's BLAS and pyarrow's jemalloc from starting threads in it, so
    each worker starts with pandas, numpy and h5py already imported,
    sharing their memory with the others as a forked worker would. Left to
    its default, the server imports the parent's main module, threads and
    all, and warns on every fork.
    """
    if not sys.platform.startswith("linux"):
        return mp.get_context("spawn")
    ctx = mp.get_context("forkserver")
    ctx.set_forkserver_preload(["pycoQC._forkserver_preload"])
    return ctx


def wait_for_workers(ps_list, error_q, poll=0.5):
    """Wait for the writer's report, raising the first error instead.

    Workers send ("error", message) when they fail and the writer sends
    ("done", report) when it has finished. A worker that dies without
    sending anything, killed for memory or failing as it starts, used to
    leave the parent waiting for ever, so the workers' exit codes are
    checked while it waits.
    """
    while True:
        try:
            kind, payload = error_q.get(timeout=poll)
        except queue.Empty:
            dead = [p for p in ps_list if p.exitcode not in (None, 0)]
            if not dead:
                continue
            # A message may have been queued just before the exit.
            try:
                kind, payload = error_q.get(timeout=poll)
            except queue.Empty:
                raise pycoQCError(
                    "A worker process stopped without reporting why (exit "
                    "code {}). If the converter is called from a script, the "
                    "call must sit under `if __name__ == \"__main__\":`, since "
                    "each worker imports the script again as it starts. "
                    "Otherwise the system may have run out of memory: try "
                    "fewer threads.".format(dead[0].exitcode))
        if kind == "done":
            return payload
        raise pycoQCError(payload)


def stop_processes(ps_list, timeout=5):
    """Make sure none of these processes outlives the call that started them.

    `Process.terminate` only sends SIGTERM; it does not wait, and it does not
    reap. A child that is blocked writing to a full queue, or one that has
    installed a SIGTERM handler of its own, can still be running when the
    parent moves on. Python then joins it from an atexit handler, so the
    parent appears to finish its work and hangs on the way out, with no
    traceback and nothing to say which child is holding it.

    coverage installs exactly such a handler when run with `sigterm = true`,
    which is how this surfaced: the test suite passed and the interpreter
    never exited.

    So terminate, wait a short while, and SIGKILL whatever is left.
    """
    for proc in ps_list:
        if proc.is_alive():
            proc.terminate()

    for proc in ps_list:
        proc.join(timeout)

    for proc in ps_list:
        if proc.is_alive():
            # Nothing is left to flush at this point, and a child that
            # ignored SIGTERM will not honour a second one either.
            proc.kill()
            proc.join(timeout)


# Held while the worker settings are in the environment, so that two
# conversions started from threads of one process cannot interleave their
# save and restore, which could leave the settings in the parent for good.
_ENV_LOCK = threading.Lock()


@contextmanager
def _worker_environment():
    """The worker settings in this process's environment while workers start.

    A spawned worker, and the fork server, inherit the environment at the
    moment they are started, and read it before they import anything. The
    parent's own numpy loaded long ago, so it is unaffected, and its values
    are put back as soon as the workers are running. A subprocess another
    thread starts in that window inherits them too.
    """
    with _ENV_LOCK:
        saved = {k: os.environ.get(k) for k in WORKER_ENV}
        os.environ.update(WORKER_ENV)
        try:
            yield
        finally:
            for key, value in saved.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def _queued_error(error_q, timeout=0.1):
    """An error left on the queue once every worker has exited, or None.

    The writer's report can reach the parent before another worker's error
    does: each process writes to a queue from a thread of its own, so the
    order of messages from two processes is not fixed. Once the workers
    have exited, everything they sent is in the pipe.
    """
    while True:
        try:
            kind, payload = error_q.get(timeout=timeout)
        except queue.Empty:
            return None
        if kind == "error":
            return payload


def _new_counts():
    return {"overall": Counter(), "fields_found": Counter(),
            "fields_not_found": Counter()}


#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class ConverterPipeline:
    """Base class for the converters, which parse their input in processes.

    A subclass calls `_configure` with the arguments every converter takes,
    stores its own, then calls `_run`. It implements:
        _work_items(): optional, yields picklable work items in the lister;
            by default the input files
        _records(item, counts): yields one dict per read, in a parser,
            optionally counting fields in counts["fields_found"] and
            counts["fields_not_found"]
        _finalise(df, notes): optional, the table just before it is written;
            strings appended to notes are logged as warnings by the parent

    The stage methods run in other processes, so they reach the converter
    as a pickled copy: anything they need must be an instance attribute set
    before `_run`, not state changed on the class in the parent.
    """

    # What the input is called in messages: "fastq", "fast5", "bam", "pod5".
    kind = None

    # Work items waiting for a parser. A large item (a block of lines) needs
    # a short queue to bound memory; a file path does not.
    in_queue_size = 1000

    # Records travel to the writer in lists of up to this many. Sending each
    # read alone spent more time pickling than the writer spent on the rest.
    RECORDS_PER_BATCH = 1000

    # Keep one row per read_id. A directory can hold a run twice, or an
    # aligned and an unaligned copy of it, and every read would be counted
    # twice; the duplicates are dropped, keeping the last file's, and the
    # parent says how many there were.
    unique_read_ids = False

    # How the error for input with no reads in it ends.
    no_reads = "hold no reads"

    # The file extensions `_input_files` looks for.
    extensions = ()

    #~~~~~~~~~~~~~~SET UP~~~~~~~~~~~~~~#
    def _configure(self, in_dir, seq_summary_fn, max_files, threads,
                   basecall_id, include_path, verbose_level):
        """Check and store the arguments every converter takes."""
        # Kept, so each worker can set it again: a worker starts with fresh
        # loggers, and printed none of its counts before this was passed on.
        self.log_level = VERBOSE_LEVELS.get(verbose_level, logging.WARNING)
        get_logger(type(self).__module__, level=self.log_level)

        self._log.info("Check input data and options")
        check_converter_io(in_dir, seq_summary_fn, threads, self.kind)

        self.in_dir = in_dir
        self.seq_summary_fn = seq_summary_fn
        self.max_files = max_files
        # One process lists the files and one writes the summary.
        self.threads = threads - 2
        self.basecall_id = basecall_id
        self.include_path = include_path
        self.verbose_level = verbose_level

        # Only a FAST5 holds numbered basecall groups to choose between. The
        # option was accepted and ignored in silence by the other converters.
        if basecall_id and self.kind != "fast5":
            self._log.warning(
                "WARNING: --basecall_id has no meaning for {} input and is "
                "ignored; it selects one of several basecall groups inside a "
                "FAST5.".format(self.kind.upper()))

    def _input_files(self):
        """Every input file under in_dir, in order, up to max_files."""
        # Counted apart from the loop: enumerate is zero based, so testing it
        # against 0 reported "no files found" for a directory holding exactly
        # one file.
        count = 0
        for fn in recursive_file_gen(dir=self.in_dir, ext=list(self.extensions)):
            yield fn
            count += 1
            if self.max_files and count >= self.max_files:
                break
        if count == 0:
            raise pycoQCError(
                "No valid {} files found in indicated folder".format(self.kind))

    #~~~~~~~~~~~~~~TO IMPLEMENT~~~~~~~~~~~~~~#
    def _work_items(self):
        return self._input_files()

    def _records(self, item, counts):
        raise NotImplementedError

    def _finalise(self, df, notes):
        return df

    def _item_file(self, item):
        """The input file a work item comes from. Items are files unless a
        converter splits them, as the FASTQ one does."""
        return item

    #~~~~~~~~~~~~~~PARENT~~~~~~~~~~~~~~#
    @property
    def _log(self):
        # Each converter logs under its own module's name, as it always has.
        return logging.getLogger(type(self).__module__)

    def _run(self):
        """Start the lister, parsers and writer, wait for them, and report."""
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

        # The writer writes the summary under its own name, so that pandas
        # infers its compression and, for a zip or gzip, the name stored in
        # it, as before; but into a directory of this run's own beside it, so
        # that nothing of the user's is touched. It takes the summary's place
        # only once every worker has finished cleanly. Set before the workers
        # start, which receive a copy of it.
        directory, name = os.path.split(os.path.abspath(self.seq_summary_fn))
        staging = tempfile.mkdtemp(prefix=".pycoqc-", dir=directory)
        self._partial_fn = os.path.join(staging, name)

        self._log.info("Start processing {} files".format(self.kind))
        try:
            with _worker_environment():
                for proc in ps_list:
                    proc.start()
            report = wait_for_workers(ps_list, error_q)
            for proc in ps_list:
                proc.join()
            late = _queued_error(error_q)
            if late is not None:
                raise pycoQCError(late)
            # Only a run every worker finished becomes the summary. The writer
            # can finish while another worker is failing, and a summary
            # missing that worker's files would pass for a complete one.
            os.replace(self._partial_fn, self.seq_summary_fn)

        except BaseException:
            self._log.info("\nAn error occurred. All processes were killed\n")
            raise

        finally:
            # Every exit from here, including a clean one, leaves nothing
            # running. Listing the exceptions to clean up after meant any
            # other kind left the whole pool alive.
            stop_processes(ps_list)
            shutil.rmtree(staging, ignore_errors=True)

        self._report(report)

    def _report(self, report):
        log = self._log
        counts = report["counts"]
        log.info("Overall counts {}".format(dict_to_str(counts["overall"])))
        if counts["fields_found"] or counts["fields_not_found"]:
            log.info("fields found {}".format(dict_to_str(counts["fields_found"])))
            log.info("fields not found {}".format(dict_to_str(counts["fields_not_found"])))
        if report["duplicates_across_files"]:
            log.warning(
                "WARNING: {:,} read id(s) appear in more than one {} file and "
                "were written once. Check whether the directory holds the same "
                "run twice.".format(report["duplicates_across_files"], self.kind))
        if report["duplicates_within_files"]:
            log.warning(
                "WARNING: {:,} read id(s) appear more than once within one {} "
                "file and were written once.".format(
                    report["duplicates_within_files"], self.kind))
        for note in report["notes"]:
            log.warning(note)
        log.warning("Total reads: {} / Average speed: {} reads/s\n".format(
            report["reads"], round(report["reads"] / report["seconds"], 2)))

    #~~~~~~~~~~~~~~WORKERS~~~~~~~~~~~~~~#
    def _start_worker(self):
        # A worker starts with fresh loggers rather than the parent's.
        return get_logger(type(self).__module__, level=self.log_level)

    def _lister(self, in_q, error_q):
        """Queue the numbered work items, then one stop signal per parser."""
        log = self._start_worker()
        try:
            n = 0
            for n, item in enumerate(self._work_items(), 1):
                in_q.put((n, item))
            log.debug("[LISTER] Queued {} work items".format(n))
        except Exception as E:
            error_q.put(("error", format_worker_error(E)))
        finally:
            for _ in range(self.threads):
                in_q.put(None)

    def _parser(self, in_q, out_q, error_q, counter_q, worker_id):
        """Turn work items into batches of records for the writer.

        Each batch carries its item's number and its own place within the
        item, so the writer can put the rows back in the order the files
        were listed, whichever parser finished first.
        """
        log = self._start_worker()
        log.debug("[PARSER_{:02}] Start".format(worker_id))
        try:
            counts = _new_counts()
            for n, item in iter(in_q.get, None):
                source = self._item_file(item)
                part = 0
                batch = []
                for record in self._records(item, counts):
                    if not record:
                        counts["overall"]["invalid reads"] += 1
                        continue
                    counts["overall"]["valid reads"] += 1
                    batch.append(record)
                    if len(batch) >= self.RECORDS_PER_BATCH:
                        out_q.put((n, part, source, batch))
                        part += 1
                        batch = []
                if batch:
                    out_q.put((n, part, source, batch))
            counter_q.put(counts)
        except Exception as E:
            error_q.put(("error", format_worker_error(E)))
        finally:
            out_q.put(None)

    def _writer(self, out_q, error_q, counter_q):
        """Collect every record, write the summary and send back a report."""
        log = self._start_worker()
        t = time()
        try:
            batches = []
            with tqdm(unit=" reads", mininterval=0.1, smoothing=0.1,
                      disable=self.verbose_level == 2) as pbar:
                for _ in range(self.threads):
                    for batch in iter(out_q.get, None):
                        batches.append(batch)
                        pbar.update(len(batch[3]))

            # Files holding no reads used to give an empty summary, which
            # pycoQC then rejected with a message about the summary rather
            # than about the input files.
            log.debug("[WRITER] Write data to file")
            batches.sort(key=lambda b: (b[0], b[1]))
            rows = [record for _, _, _, batch in batches for record in batch]
            if not rows:
                raise pycoQCError(
                    "No reads found in the {} file(s). The files were read "
                    "but {}.".format(self.kind, self.no_reads))
            # The file each row came from, which a split file shares between
            # several work items.
            sources = [source for _, _, source, batch in batches for _ in batch]
            df = pd.DataFrame(rows)
            del rows, batches

            across = within = 0
            if self.unique_read_ids and "read_id" in df:
                repeated = df.duplicated("read_id", keep="last")
                if repeated.any():
                    # A read twice in one file is a different problem from
                    # a run present twice in the directory, and only the
                    # second is worth pointing the user at their files for.
                    within = int(pd.DataFrame(
                        {"read_id": df["read_id"].to_numpy(), "source": sources}
                    ).duplicated().sum())
                    across = int(repeated.sum()) - within
                    df = df[~repeated].reset_index(drop=True)

            notes = []
            df = self._finalise(df, notes)
            df.to_csv(self._partial_fn, sep="\t", index=False)

            counts = _new_counts()
            for _ in range(self.threads):
                for key, counter in counter_q.get().items():
                    counts[key].update(counter)

            error_q.put(("done", {
                "reads": len(df), "seconds": time() - t, "counts": counts,
                "duplicates_across_files": across,
                "duplicates_within_files": within, "notes": notes}))

        except Exception as E:
            error_q.put(("error", format_worker_error(E)))

# -*- coding: utf-8 -*-
"""The process pipeline every converter shares.

A lister process queues work items, a pool of parser processes turns each
item into one record per read, and a writer process collects the records
into the summary file:

    lister --in_q--> parser x N --out_q--> writer
      |                  |                     |
      '----- a status pipe each, to the parent ------'

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
import shutil
import sys
import tempfile
import threading
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from multiprocessing.connection import wait as wait_ready
from time import time

# Third party imports
import numpy as np
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


# Where an error came from, earliest stage first. When several workers fail,
# the earliest stage's error is the cause and the others its consequences: a
# lister that finds no files leaves the writer with no reads to write.
LISTER, PARSER, WRITER, SILENT = 0, 1, 2, 3


_SILENT_EXIT = (
    "A worker process stopped without reporting why (exit code {}). If the "
    "converter is called from a script, the call must sit under "
    "`if __name__ == \"__main__\":`, since each worker imports the script again "
    "as it starts. Otherwise the system may have run out of memory: try fewer "
    "threads.")


# The parent waits on two handles per worker, its pipe and its exit. Before
# Python 3.13, Windows waits on at most 63 at once, and a run of 32 or more
# workers failed as soon as it started; 3.13 waits on more in batches.
WINDOWS_MAX_THREADS = 31


def usable_threads(threads, platform=sys.platform, version=sys.version_info):
    """The threads a conversion can use: `threads`, or fewer where this
    Python cannot wait on that many workers. Raises below the minimum.

    One process lists the files and one writes the summary, so at least
    three are needed for any parsing to happen.
    """
    if threads < 3:
        raise pycoQCError("At least 3 threads required")
    if platform == "win32" and version < (3, 13):
        return min(threads, WINDOWS_MAX_THREADS)
    return threads


def name_start(name, limit=50):
    """The start of a file name, at most `limit` bytes, whole characters only.

    A staging directory is named after the summary so that one left by a
    killed run is easy to find, and only the start of the name is needed
    for that. The limit is in bytes: a file name holds at most 255 of them,
    and 150 characters of Chinese is 450. Bytes that are not valid in the
    file system's encoding are kept, as one byte each.

    The summary goes inside the directory under its whole name, so the path
    written to is up to about 72 characters longer than the summary's own.
    A short start keeps that small, but on Windows without long paths
    enabled, a summary within that much of the 260 character limit cannot
    be written.
    """
    start, size = [], 0
    for char in name:
        size += len(os.fsencode(char))
        if size > limit:
            break
        start.append(char)
    return "".join(start)


def wait_for_workers(workers, grace=0.5):
    """Wait for every worker to finish, then return the writer's report or
    raise the error that explains the failure.

    * workers
        (process, connection) pairs. Each worker reports on a pipe of its
        own: ("error", stage, message) when it fails, and the writer
        ("done", report) when it has finished.
    * grace
        Seconds the others are given, once an error has arrived, before
        the error is raised: a worker downstream of a failure may never
        exit on its own.

    The parent sleeps until a worker sends something or exits, and reads
    each pipe to its end, so a worker is never left unable to exit with a
    message unread. When a worker has exited, everything it sent is read
    before its exit code is looked at.

    Of several errors, the earliest stage's is raised (lister, then
    parsers, then writer), since the later ones are usually its
    consequences: a lister that finds no files leaves the writer with no
    reads. A worker that dies without reporting, killed for memory or
    failing as it starts, is noticed from its exit code. Once the writer
    has reported, though, every file has been read and the summary
    written, so a worker crashing on its way out does not fail the run.
    """
    report = None
    errors = []
    deadline = None
    pipes = {conn for _, conn in workers}
    sentinels = {proc.sentinel: (proc, conn) for proc, conn in workers}

    def take(message):
        nonlocal report, deadline
        if message[0] == "done":
            report = message[1]
            return
        errors.append(message[1:])
        if deadline is None:
            deadline = time() + grace

    def read_to_end(conn):
        # Only called on a worker that has exited: its end is closed, so
        # this returns at the end of what it sent rather than waiting.
        while conn in pipes:
            try:
                take(conn.recv())
            except EOFError:
                pipes.discard(conn)

    while pipes or sentinels:
        timeout = None
        if errors:
            timeout = deadline - time()
            if timeout <= 0:
                break
        # Each object once: on Windows, CPython 3.13's batched wait can list
        # a handle twice, and a second recv on a drained pipe blocks.
        ready = list(dict.fromkeys(wait_ready(list(pipes) + list(sentinels), timeout)))
        # Messages first, so a report is in hand before an exit is judged.
        for obj in sorted(ready, key=lambda o: o in sentinels):
            if obj in pipes:
                try:
                    take(obj.recv())
                except EOFError:
                    pipes.discard(obj)
            elif obj in sentinels:
                proc, conn = sentinels.pop(obj)
                proc.join()
                read_to_end(conn)
                if proc.exitcode != 0 and report is None:
                    take(("error", SILENT, _SILENT_EXIT.format(proc.exitcode)))

    if errors:
        # min keeps the first to arrive among errors of the same stage.
        raise pycoQCError(min(errors, key=lambda e: e[0])[1])
    if report is None:
        raise pycoQCError("The workers finished without writing a summary")
    return report


def _send_error(status, stage, error):
    """Report a worker's failure. An Exception reported is done with;
    anything else, or an Exception that could not be reported, is raised
    again, so the worker exits non-zero.

    Only an Exception was caught. A lister stopped by SystemExit or
    KeyboardInterrupt still sent the parsers their stop signals on its way
    out, so the writer wrote and reported a summary missing the files not yet
    listed. Its exit code came after the report, which is not looked at, and
    the short summary was kept as complete.
    """
    try:
        status.send(("error", stage, format_worker_error(error)))
        sent = True
    except OSError:
        sent = False
    # A send fails only once the parent has stopped listening: on Ctrl-C, or
    # when it is already raising another worker's error. The parent is past
    # waiting by then and will not keep a summary, so this is about what the
    # worker leaves behind: its own error, not the pipe's, which would hide
    # it, and a non-zero exit, not a clean one. Raised here rather than in
    # the except clause, so that it does not carry the pipe's error, and the
    # pickled message, along.
    if not sent or not isinstance(error, Exception):
        raise error


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
    # A process that was never started has nothing to stop, and joining it
    # raises: start() itself can fail partway through the list.
    ps_list = [proc for proc in ps_list if proc.pid is not None]

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
        _item_file(item): the input file a work item comes from; needed by a
            converter whose items are not plain file paths, so that read ids
            repeated within one file are not reported as across files

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

    # The fields written when none are named: fields=None or an empty list in
    # Python, and the --fields default of the command line, which reads it
    # from here, as the notebook help does. A converter, or a subclass of
    # one, sets its default here and nowhere else; its __init__ takes
    # fields=None. None here is every field the input carries. The command
    # line once kept lists of its own, and calling the FASTQ class directly
    # silently dropped the pass or fail column the command produced.
    DEFAULT_FIELDS = None

    @classmethod
    def help_defaults(cls):
        """Defaults to show in place of the signature's: fields=None stands
        for DEFAULT_FIELDS, listed, or every field where that is None. The
        command line and the notebook help both read them from here."""
        return {"fields": None if cls.DEFAULT_FIELDS is None else list(cls.DEFAULT_FIELDS)}

    #~~~~~~~~~~~~~~SET UP~~~~~~~~~~~~~~#
    def _configure(self, in_dir, seq_summary_fn, max_files, threads,
                   basecall_id, include_path, verbose_level):
        """Check and store the arguments every converter takes."""
        # Kept, so each worker can set it again: a worker starts with fresh
        # loggers, and printed none of its counts before this was passed on.
        self.log_level = VERBOSE_LEVELS.get(verbose_level, logging.WARNING)
        get_logger(type(self).__module__, level=self.log_level)

        self._log.info("Check input data and options")
        check_converter_io(in_dir, seq_summary_fn, self.kind)
        usable = usable_threads(threads)
        if usable < threads:
            self._log.warning(
                "WARNING: --threads lowered from {} to {}: Python {}.{} on "
                "Windows waits on at most 63 handles, two per thread. Python "
                "3.13 and later have no such limit.".format(
                    threads, usable, sys.version_info[0], sys.version_info[1]))

        self.in_dir = in_dir
        self.seq_summary_fn = seq_summary_fn
        self.max_files = max_files
        # One process lists the files and one writes the summary.
        self.threads = usable - 2
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

    def _requested_fields(self, fields, keep_path=False):
        """--fields, each named once, and without path unless kept.

        A field named twice was read and counted twice, and written as two
        columns by the BAM converter.

        The path column is the file a read came from, written with
        --include_path. Asked for without it, each converter blamed the input:
        the BAM one said a BAM carries no such field, the FASTQ one counted it
        missing from every read, and the FAST5 one rejected it as invalid.
        With --include_path, `keep_path` leaves it where it was listed, for a
        converter that places the columns in the order asked for.

        None or an empty list is the converter's DEFAULT_FIELDS, as both
        were for BAM: the FASTQ converter wrote nothing for each read of an
        empty list and failed on finding no reads, and with FAST5 failed on
        None. Asking for path alone, without --include_path, leaves the read
        ids, not every field and not none.

        See _field_list for what is accepted as names.
        """
        fields = self._field_list(fields)
        if not fields:
            fields = self.DEFAULT_FIELDS
        if fields is None:
            return None
        fields = list(dict.fromkeys(fields))
        if "path" in fields:
            if not self.include_path:
                self._log.warning(
                    "WARNING: the path column is only written with "
                    "--include_path, so path is left out of the summary")
            if not (keep_path and self.include_path):
                fields.remove("path")
        return fields or ["read_id"]

    @staticmethod
    def _field_list(fields):
        """`fields` as a list of plain, stripped names, or None.

        Accepted: one name as a string, where a blank one is no name at all,
        or names in a list, a tuple, a one dimensional numpy array, a pandas
        Index or an iterator such as a generator. Each name must be a string
        with something in it. These are named rather than anything iterable
        being taken: a string was split one field per letter, bytes one per
        byte, a set gave its columns in a different order on each run, and a
        dict or a Series gave its keys or its values without a word.

        An error raised by the caller's own iterator is left as it is.
        """
        if fields is None:
            return None
        if isinstance(fields, str):
            names = [fields] if fields.strip() else []
        elif (isinstance(fields, (list, tuple, pd.Index, Iterator))
              or (isinstance(fields, np.ndarray) and fields.ndim == 1)):
            names = list(fields)
        else:
            raise pycoQCError(
                "fields must be a field name, or field names in a list, a "
                "tuple, an array, an Index or an iterator, not {}".format(
                    type(fields).__name__))
        for name in names:
            if not isinstance(name, str) or not name.strip():
                raise pycoQCError(
                    "fields must be field names, as strings with something "
                    "in them; got {!r}".format(name))
        # numpy's str_ is a str, but plain names keep the columns plain.
        return [str(name).strip() for name in names]

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
        counter_q = ctx.Queue()

        # Each worker reports on a pipe of its own, which the parent can
        # sleep on together with the workers' exits. A worker is paired with
        # its pipe as it is made, and given its end as its first argument.
        workers = []

        def worker(target, *args):
            reader, writer = ctx.Pipe(duplex=False)
            workers.append((ctx.Process(target=target, args=(writer,) + args),
                            reader, writer))

        worker(self._lister, in_q)
        for i in range(self.threads):
            worker(self._parser, in_q, out_q, counter_q, i)
        worker(self._writer, out_q, counter_q)
        ps_list = [proc for proc, _, _ in workers]

        self._log.info("Start processing {} files".format(self.kind))
        staging = None
        try:
            # The writer writes the summary under its own name, so that
            # pandas infers its compression and, for a zip or gzip, the name
            # stored in it, as before; but into a directory of this run's
            # own beside it, so that nothing of the user's is touched. It
            # takes the summary's place only once every worker has finished
            # cleanly. The directory is named after the summary and not
            # hidden, so that one left by a killed run is easy to find. Set
            # before the workers start, which receive a copy of the path.
            directory, name = os.path.split(os.path.abspath(self.seq_summary_fn))
            # Only the start of the name: see name_start.
            staging = tempfile.mkdtemp(prefix=name_start(name) + ".pycoqc-tmp-",
                                       dir=directory)
            self._partial_fn = os.path.join(staging, name)

            with _worker_environment():
                for proc in ps_list:
                    proc.start()
            # The workers hold their own ends now. Closing the parent's
            # copies is what lets a pipe report its end once its worker has
            # exited.
            for _, _, writer in workers:
                writer.close()
            report = wait_for_workers(
                [(proc, reader) for proc, reader, _ in workers])
            for proc in ps_list:
                proc.join()
            # Only a run every worker finished becomes the summary. The writer
            # can finish while another worker is failing, and a summary
            # missing that worker's files would pass for a complete one.
            os.replace(self._partial_fn, self.seq_summary_fn)

        except BaseException:
            # Only when there were any: creating the temporary directory
            # fails before a single worker has started.
            if any(proc.pid is not None for proc in ps_list):
                self._log.info("\nAn error occurred. All processes were killed\n")
            raise

        finally:
            # Every exit from here, including a clean one, leaves nothing
            # running. Listing the exceptions to clean up after meant any
            # other kind left the whole pool alive.
            stop_processes(ps_list)
            for _, reader, writer in workers:
                reader.close()
                writer.close()
            if staging is not None:
                shutil.rmtree(staging, ignore_errors=True)
                if os.path.exists(staging):
                    self._log.warning(
                        "WARNING: could not remove the temporary directory "
                        "{}".format(staging))

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

    def _lister(self, status, in_q):
        """Queue the numbered work items, then one stop signal per parser."""
        log = self._start_worker()
        try:
            n = 0
            for n, item in enumerate(self._work_items(), 1):
                in_q.put((n, item))
            log.debug("[LISTER] Queued {} work items".format(n))
        except BaseException as E:
            _send_error(status, LISTER, E)
        finally:
            for _ in range(self.threads):
                in_q.put(None)

    def _parser(self, status, in_q, out_q, counter_q, worker_id):
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
                # Only needed to tell repeated read ids apart by file.
                source = self._item_file(item) if self.unique_read_ids else None
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
        except BaseException as E:
            _send_error(status, PARSER, E)
        finally:
            out_q.put(None)

    def _writer(self, status, out_q, counter_q):
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
            # several work items. Only kept where repeats are looked for.
            sources = ([source for _, _, source, batch in batches for _ in batch]
                       if self.unique_read_ids else None)
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

            status.send(("done", {
                "reads": len(df), "seconds": time() - t, "counts": counts,
                "duplicates_across_files": across,
                "duplicates_within_files": within, "notes": notes}))

        except BaseException as E:
            _send_error(status, WRITER, E)

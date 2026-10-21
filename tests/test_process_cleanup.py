# -*- coding: utf-8 -*-
"""The converters must not leave worker processes behind.

`Process.terminate` only sends SIGTERM, and the error path used to send it
without ever waiting for the children to go. Anything that survived the
signal was left for Python's atexit handler to join, so the work looked
finished and the interpreter hung on the way out instead, with no traceback
naming the child that held it. CI caught it as a test run that reported
"427 passed" and then sat there until it was cancelled 20 minutes later,
with two orphaned python processes in the cleanup log.

A child can survive SIGTERM by handling it, which is exactly what coverage
does when configured with `sigterm = true`.
"""

import multiprocessing as mp

import pytest

from pycoQC.common import pycoQCError, stop_processes
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary


def _survivors():
    """Children still running, having reaped the ones that have exited."""
    return [p for p in mp.active_children() if p.is_alive()]


class TestNothingIsLeftRunning:

    def test_fast5_error_path_leaves_no_children(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(pycoQCError):
            Fast5_to_seq_summary(fast5_dir=str(empty),
                                 seq_summary_fn=str(tmp_path / "x.txt"),
                                 threads=3)
        assert _survivors() == []

    def test_fastq_error_path_leaves_no_children(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(pycoQCError):
            Fastq_to_seq_summary(fastq_dir=str(empty),
                                 seq_summary_fn=str(tmp_path / "x.txt"),
                                 threads=3)
        assert _survivors() == []

    def test_fast5_success_path_leaves_no_children(
            self, tmp_path, basecalled_fast5_dir):
        Fast5_to_seq_summary(fast5_dir=str(basecalled_fast5_dir),
                             seq_summary_fn=str(tmp_path / "s.txt"),
                             threads=3)
        assert _survivors() == []

    def test_fastq_success_path_leaves_no_children(self, tmp_path, fastq_dir):
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir),
                             seq_summary_fn=str(tmp_path / "s.txt"),
                             threads=3)
        assert _survivors() == []


class TestStopProcesses:

    def test_a_child_ignoring_sigterm_is_still_stopped(self):
        """coverage's SIGTERM handler is the real case this stands in for."""
        ready = mp.Event()
        ps = mp.Process(target=_ignores_sigterm, args=(ready,))
        ps.start()
        try:
            # Without this the child is usually killed before it installs
            # the handler, and the test passes whatever the cleanup does.
            assert ready.wait(30), "child never got going"
            stop_processes([ps], timeout=2)
            assert not ps.is_alive()
        finally:
            if ps.is_alive():       # the assertion above failed
                ps.kill()
                ps.join(5)

    def test_an_already_finished_child_is_not_signalled(self):
        ps = mp.Process(target=_returns_immediately)
        ps.start()
        ps.join(10)
        stop_processes([ps], timeout=2)
        assert ps.exitcode == 0     # not -SIGTERM, and not -SIGKILL


def _ignores_sigterm(ready):
    import signal
    import time
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    ready.set()
    time.sleep(120)


def _returns_immediately():
    return

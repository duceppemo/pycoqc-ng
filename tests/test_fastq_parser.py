# -*- coding: utf-8 -*-
"""The chunked fastq reader.

The reader hands blocks of lines to worker processes, so a block that does
not land on a read boundary, or a thread count that changes how blocks are
distributed, must not alter the result.
"""

import gzip

import pandas as pd
import pytest

from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary, FastqParser
from conftest import N_READS, reference_mean_qscore


class TestHelpers:

    @pytest.mark.parametrize("path, expected", [
        ("/data/fastq_pass/reads.fastq.gz", "pass"),
        ("/data/fastq_fail/reads.fastq.gz", "fail"),
        ("/data/reads.fastq", "pass"),
        ("/data/fail/reads.fastq", "fail"),
        ("/data/run_fail_01/reads.fastq", "fail"),
        # The basecaller writes the word "fail" on its own, between
        # separators. A directory that merely starts with those letters is
        # not a statement about the reads, and treating it as one marked
        # entire good runs as failing.
        ("/mnt/failsafe_backup/fastq_pass/reads.fastq.gz", "pass"),
        ("/mnt/failover/project/fastq_pass/reads.fastq.gz", "pass"),
        ("/data/failed_run/reads.fastq", "pass"),
    ])
    def test_pass_fail_flag_from_path(self, path, expected):
        assert FastqParser.pass_fail_flag(path) == expected

    def test_open_fastq_handles_gzip_and_plain(self, tmp_path):
        plain = tmp_path / "a.fastq"
        plain.write_text("@r1\nACGT\n+\nIIII\n")
        gzipped = tmp_path / "b.fastq.gz"
        with gzip.open(gzipped, "wt") as fh:
            fh.write("@r1\nACGT\n+\nIIII\n")

        for fn in (plain, gzipped):
            with FastqParser.open_fastq(str(fn)) as fh:
                assert fh.readline().startswith("@r1")

    def test_make_chunks_covers_every_line_exactly_once(self):
        lines = [f"line{i}\n" for i in range(1000)]
        for size in (1, 3, 4, 7, 400, 4000):
            chunks = list(FastqParser.make_chunks(iter(lines), size))
            assert [x for chunk in chunks for x in chunk] == lines
            assert all(len(c) <= size for c in chunks)

    def test_read_entry_parses_a_whole_block(self, reads):
        block = []
        for read_id, seq, qual in reads[:5]:
            block += [f"@{read_id} runid=r1 ch=1 start_time=2022-06-27T23:00:28Z\n",
                      seq + "\n", "+\n", qual + "\n"]
        parsed = FastqParser.read_entry(block, "pass")
        assert len(parsed) == 5
        for read_id, seq, qual in reads[:5]:
            assert parsed[read_id]["sequence_length_template"] == len(seq)
            assert parsed[read_id]["mean_qscore_template"] == pytest.approx(
                round(reference_mean_qscore(qual), 2), abs=0.011)


class TestChunking:

    @pytest.fixture
    def summary_for(self, tmp_path, fastq_dir):
        def build(threads, reads_per_chunk=None, tag=""):
            original = Fastq_to_seq_summary.READS_PER_CHUNK
            if reads_per_chunk is not None:
                Fastq_to_seq_summary.READS_PER_CHUNK = reads_per_chunk
            try:
                out = tmp_path / f"summary_{threads}_{reads_per_chunk}{tag}.txt"
                Fastq_to_seq_summary(fastq_dir=str(fastq_dir),
                                     seq_summary_fn=str(out), threads=threads)
                return pd.read_csv(out, sep="\t").sort_values("read_id")
            finally:
                Fastq_to_seq_summary.READS_PER_CHUNK = original
        return build

    @pytest.mark.parametrize("reads_per_chunk", [1, 2, 7, 101, 1000])
    def test_all_reads_survive_any_chunk_size(self, summary_for, reads_per_chunk):
        """Chunk size is counted in reads, so no block can split a record no
        matter how the size is set."""
        df = summary_for(threads=3, reads_per_chunk=reads_per_chunk)
        assert len(df) == N_READS
        assert df["read_id"].nunique() == N_READS

    def test_result_is_independent_of_thread_count(self, summary_for):
        """Blocks are distributed differently as the worker count changes, so
        this pins that the output does not depend on it."""
        few = summary_for(threads=3, tag="a")
        many = summary_for(threads=8, tag="b")
        pd.testing.assert_frame_equal(
            few.reset_index(drop=True), many.reset_index(drop=True))

    def test_values_correct_when_split_across_many_chunks(
            self, summary_for, reads):
        df = summary_for(threads=4, reads_per_chunk=10).set_index("read_id")
        for read_id, seq, qual in reads:
            assert df.loc[read_id, "sequence_length_template"] == len(seq)
            assert df.loc[read_id, "mean_qscore_template"] == pytest.approx(
                round(reference_mean_qscore(qual), 2), abs=0.011)


class TestTruncatedInput:

    def test_a_block_ending_mid_record_warns_and_skips(self, caplog):
        """Reads are only ever cut short by a truncated file, which must be
        reported rather than quietly losing the record."""
        block = ["@r1 runid=x ch=1 start_time=2022-06-27T23:00:28Z\n",
                 "ACGT\n", "+\n", "IIII\n",
                 "@r2 runid=x ch=1 start_time=2022-06-27T23:00:28Z\n",
                 "ACGT\n"]
        with caplog.at_level("WARNING"):
            parsed = FastqParser.read_entry(block, "pass")
        assert set(parsed) == {"r1"}
        assert "Truncated fastq record" in caplog.text

    def test_complete_blocks_do_not_warn(self, caplog):
        block = ["@r1 runid=x ch=1 start_time=2022-06-27T23:00:28Z\n",
                 "ACGT\n", "+\n", "IIII\n"]
        with caplog.at_level("WARNING"):
            parsed = FastqParser.read_entry(block, "pass")
        assert set(parsed) == {"r1"}
        assert "Truncated" not in caplog.text

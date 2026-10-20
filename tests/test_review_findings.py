# -*- coding: utf-8 -*-
"""Defects found in review, each pinned so it cannot come back.

Most of these produced plausible but wrong output rather than failing, which
is the dangerous kind for a tool whose numbers end up in a methods section.
"""

import gzip
import random

import pandas as pd
import pytest

from pycoQC.common import check_arg
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot


def write_summary(path, n=2000, run_id="run1", barcodes=True, seed=3,
                  **extra):
    rng = random.Random(seed)
    df = pd.DataFrame({
        "read_id": [f"r{i:06d}" for i in range(n)],
        "run_id": run_id,
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
        "sequence_length_template": [rng.randint(200, 5000) for _ in range(n)],
        "mean_qscore_template": [round(rng.uniform(4, 14), 2) for _ in range(n)],
    })
    if barcodes:
        df["barcode_arrangement"] = [f"barcode{i % 4 + 1:02d}" for i in range(n)]
    for k, v in extra.items():
        df[k] = v
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture(scope="module")
def plotter(tmp_path_factory):
    path = write_summary(tmp_path_factory.mktemp("bc") / "s.txt", n=2000)
    parser = pycoQC_parse(summary_file=str(path), quiet=True)
    # A tenth of the reads, so an unscaled count is obviously wrong.
    return parser, pycoQC_plot(parser, sample=200, quiet=True)


class TestBarcodeCountsAreScaledToTheRun:
    """Plots are built from a sample of at most `sample` reads. Counts taken
    from it must be scaled back up, as the over-time plots already did."""

    def test_per_barcode_reads_sum_to_the_run(self, plotter):
        parser, p = plotter
        cells = p.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        reads = [v for v in cells[1] if isinstance(v, (int, float))]
        assert sum(reads) == pytest.approx(len(parser.reads_df), rel=0.02)

    def test_per_barcode_bases_are_scaled_too(self, plotter):
        parser, p = plotter
        cells = p.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        bases = [v for v in cells[2] if isinstance(v, (int, float))]
        assert sum(bases) == pytest.approx(
            int(parser.reads_df["read_len"].sum()), rel=0.05)

    def test_pie_counts_sum_to_the_run(self, plotter):
        parser, p = plotter
        for trace in p.barcode_counts()["data"]:
            values = getattr(trace, "values", None)
            if values is not None:
                assert sum(values) == pytest.approx(
                    len(parser.reads_df), rel=0.02)
                return
        pytest.fail("no pie trace found")

    def test_unsampled_run_is_unchanged(self, tmp_path):
        """With no sampling the factor is 1 and counts must be exact."""
        path = write_summary(tmp_path / "s.txt", n=600)
        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        p = pycoQC_plot(parser, sample=0, quiet=True)
        cells = p.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        reads = [v for v in cells[1] if isinstance(v, (int, float))]
        assert sum(reads) == 600


class TestNumberFormats:

    def test_alignment_summary_formats(self, aligned_bam, tmp_path):
        """d3's "r" rounds to significant digits, so ".1r" turned a median
        read length of 4,732 into 5,000, and identity is a 0-1 fraction that
        ".1f" flattened to 1.0."""
        path = write_summary(tmp_path / "s.txt", n=600, barcodes=False)
        df = pd.read_csv(path, sep="\t")
        df["read_id"] = [f"r{i:05d}" for i in range(len(df))]
        df.to_csv(path, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(path),
                              bam_file=str(aligned_bam), quiet=True)
        fig = pycoQC_plot(parser, quiet=True).alignment_summary()
        fmt = fig["data"][0]["cells"]["format"]
        assert fmt[5] == ","      # median read length, a count of bases
        assert fmt[6] == ".4f"    # identity, a fraction between 0 and 1


class TestBamConverterRobustness:

    def _bam(self, path, records):
        pysam = pytest.importorskip("pysam")
        header = {"HD": {"VN": "1.6", "SO": "unknown"},
                  "SQ": [{"SN": "chr1", "LN": 10000}],
                  "RG": [{"ID": "rg1", "DS": "runid=r1 basecall_model=m1",
                          "PU": "FC1", "PM": "minion", "DT": "2026-01-01T00:00:00Z"}]}
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            for rec in records(pysam):
                bam.write(rec)
        return path

    def test_secondary_and_supplementary_are_skipped(self, tmp_path):
        """Only the primary record carries the whole read; a supplementary
        one is hard clipped and overwrote the read with a fragment."""
        from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary

        def records(pysam):
            primary = pysam.AlignedSegment()
            primary.query_name = "read1"
            primary.query_sequence = "ACGT" * 25          # 100 bases
            primary.query_qualities = pysam.qualitystring_to_array("I" * 100)
            primary.flag = 0
            primary.reference_id = 0
            primary.reference_start = 10
            primary.cigartuples = [(0, 100)]
            primary.set_tags([("RG", "rg1", "Z"), ("ch", 1, "i"),
                              ("st", "2026-01-01T00:00:00Z", "Z")])

            supp = pysam.AlignedSegment()
            supp.query_name = "read1"
            supp.query_sequence = "ACGT" * 10             # clipped fragment
            supp.query_qualities = pysam.qualitystring_to_array("I" * 40)
            supp.flag = 2048
            supp.reference_id = 0
            supp.reference_start = 500
            supp.cigartuples = [(0, 40), (5, 60)]
            supp.set_tags([("RG", "rg1", "Z"), ("ch", 1, "i"),
                           ("st", "2026-01-01T00:00:00Z", "Z")])
            return [primary, supp]

        d = tmp_path / "bam"
        d.mkdir()
        self._bam(d / "r.bam", records)
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 1
        assert df["sequence_length_template"].iloc[0] == 100

    def test_record_without_a_sequence_does_not_crash(self, tmp_path):
        """A secondary record often carries no SEQ at all."""
        from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary

        def records(pysam):
            primary = pysam.AlignedSegment()
            primary.query_name = "read1"
            primary.query_sequence = "ACGT" * 25
            primary.query_qualities = pysam.qualitystring_to_array("I" * 100)
            primary.flag = 0
            primary.reference_id = 0
            primary.reference_start = 10
            primary.cigartuples = [(0, 100)]
            primary.set_tags([("RG", "rg1", "Z"), ("ch", 1, "i"),
                              ("st", "2026-01-01T00:00:00Z", "Z")])

            secondary = pysam.AlignedSegment()
            secondary.query_name = "read1"
            secondary.flag = 256
            secondary.reference_id = 0
            secondary.reference_start = 900
            secondary.cigartuples = [(0, 100)]
            secondary.set_tags([("RG", "rg1", "Z"), ("ch", 1, "i"),
                                ("st", "2026-01-01T00:00:00Z", "Z")])
            return [primary, secondary]

        d = tmp_path / "bam"
        d.mkdir()
        self._bam(d / "r.bam", records)
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(out), threads=3)
        assert len(pd.read_csv(out, sep="\t")) == 1

    def test_read_without_a_matching_read_group(self, tmp_path):
        """The blanks the docstring promises, not a KeyError."""
        from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
        pysam = pytest.importorskip("pysam")

        d = tmp_path / "bam"
        d.mkdir()
        header = {"HD": {"VN": "1.6", "SO": "unknown"}}
        with pysam.AlignmentFile(str(d / "r.bam"), "wb", header=header) as bam:
            rec = pysam.AlignedSegment()
            rec.query_name = "read1"
            rec.query_sequence = "ACGT" * 25
            rec.query_qualities = pysam.qualitystring_to_array("I" * 100)
            rec.flag = 4
            rec.set_tags([("ch", 1, "i"), ("st", "2026-01-01T00:00:00Z", "Z")])
            bam.write(rec)

        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 1


class TestFastqConverterRobustness:

    def test_zero_length_read_does_not_abort_the_run(self, tmp_path):
        """An empty record is valid fastq, and dividing by its length killed
        the conversion of every remaining file."""
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        d = tmp_path / "fq"
        d.mkdir()
        with gzip.open(d / "r.fastq.gz", "wt") as fh:
            fh.write("@empty runid=r1 ch=1 start_time=2026-01-01T00:00:00Z\n\n+\n\n")
            for i in range(5):
                fh.write(f"@read{i} runid=r1 ch=1 start_time=2026-01-01T00:0{i}:00Z\n"
                         "ACGTACGTAC\n+\nIIIIIIIIII\n")
        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(fastq_dir=str(d), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 6
        assert (df["gc_percent"] >= 0).all()

    def test_fields_without_start_time(self, tmp_path):
        """--fields lets a caller drop start_time, which used to raise
        KeyError after every file had already been parsed."""
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        d = tmp_path / "fq"
        d.mkdir()
        with gzip.open(d / "r.fastq.gz", "wt") as fh:
            for i in range(5):
                fh.write(f"@read{i} runid=r1 ch=1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(
            fastq_dir=str(d), seq_summary_fn=str(out), threads=3,
            fields=["read_id", "run_id", "sequence_length_template",
                    "mean_qscore_template"])
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 5
        assert "start_time" not in df.columns


class TestParserRobustness:

    def test_summary_with_barcodes_plus_a_barcode_file(self, tmp_path):
        """Guppy writes barcode_arrangement into the summary when barcoding
        is on; merging a separate file then left barcode_x and barcode_y."""
        summary = write_summary(tmp_path / "s.txt", n=600, barcodes=True)
        barcodes = tmp_path / "bc.txt"
        pd.DataFrame({
            "read_id": [f"r{i:06d}" for i in range(600)],
            "barcode_arrangement": ["barcode09"] * 600,
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        # The separate file is the more specific source and wins.
        assert set(parser.reads_df["barcode"]) == {"barcode09"}

    def test_unknown_runid_does_not_poison_later_runs(self, tmp_path):
        """max() of an empty selection is NaN, and the offset carried it into
        every later run's start_time without raising."""
        first = write_summary(tmp_path / "a.txt", n=400, run_id="RUN_A")
        second = write_summary(tmp_path / "b.txt", n=400, run_id="RUN_B",
                               seed=9)
        parser = pycoQC_parse(
            summary_file=[str(first), str(second)],
            runid_list=["RUN_A", "TYPO", "RUN_B"], quiet=True)
        assert parser.reads_df["start_time"].notna().all()
        assert set(parser.reads_df["run_id"]) == {"RUN_A", "RUN_B"}

    def test_1d2_summary_with_both_length_columns(self, tmp_path):
        """A 1D2 file carries the template and the 2d column, and renaming
        both onto read_len left duplicate labels pandas refuses to reindex."""
        path = write_summary(tmp_path / "s.txt", n=600, barcodes=False)
        df = pd.read_csv(path, sep="\t")
        df["sequence_length_2d"] = df["sequence_length_template"] * 2
        df["mean_qscore_2d"] = df["mean_qscore_template"]
        df.to_csv(path, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        assert len(parser.reads_df) == 600
        assert list(parser.reads_df.columns).count("read_len") == 1


class TestCheckArgBounds:

    def test_a_minimum_of_zero_is_applied(self):
        """0 is falsy, so the bound was skipped and negatives were accepted."""
        with pytest.raises(Exception, match="too low"):
            check_arg("min_pass_len", -5, int, min=0)

    def test_a_maximum_of_zero_is_applied(self):
        with pytest.raises(Exception, match="too high"):
            check_arg("n", 5, int, max=0)

    def test_zero_itself_still_passes(self):
        assert check_arg("n", 0, int, min=0) == 0

    def test_negative_min_pass_len_is_rejected(self, tmp_path):
        from pycoQC.pycoQC import pycoQC
        summary = write_summary(tmp_path / "s.txt", n=600)
        with pytest.raises(Exception):
            pycoQC(summary_file=str(summary),
                   html_outfile=str(tmp_path / "r.html"),
                   min_pass_len=-1, quiet=True)


class TestConfigWithoutPlotTitle:

    def test_entry_without_a_title_is_tolerated(self, tmp_path, guppy_summary):
        """A hand written config entry need not carry plot_title."""
        import json
        from pycoQC.pycoQC import pycoQC
        config = tmp_path / "c.json"
        config.write_text(json.dumps({"read_len_1D": {}}))
        html = tmp_path / "r.html"
        pycoQC(summary_file=str(guppy_summary), html_outfile=str(html),
               config_file=str(config), quiet=True)
        assert html.exists()


@pytest.fixture(scope="module")
def coverage_fig(tmp_path_factory, aligned_bam):
    path = write_summary(tmp_path_factory.mktemp("cov") / "s.txt",
                         n=600, barcodes=False)
    df = pd.read_csv(path, sep="\t")
    df["read_id"] = [f"r{i:05d}" for i in range(len(df))]
    df.to_csv(path, sep="\t", index=False)
    parser = pycoQC_parse(summary_file=str(path),
                          bam_file=str(aligned_bam), quiet=True)
    return pycoQC_plot(parser, quiet=True).alignment_coverage(smooth_sigma=0)


class TestPlotBinningAlignment:

    def test_coverage_x_and_y_are_the_same_length(self, coverage_fig):
        """bincount returned however many bins happened to be populated, so
        the trace could be shorter than its axis."""
        trace = coverage_fig["data"][0]
        assert len(trace["x"]) == len(trace["y"])

    def test_first_bin_can_hold_coverage(self, coverage_fig):
        """digitize against bins starting at zero never returned 0, so the
        first bin was always empty and every bin drew one place right of the
        reference it belongs to."""
        import numpy as np
        y = np.asarray(coverage_fig["data"][0]["y"], dtype=float)
        assert y[0] > 0

    def test_channel_heatmap_axis_matches_its_rows(self, guppy_summary):
        """y was bins[1:], one shorter than z, and plotly then treats the
        values as brick centres and draws every row one bin late."""
        import numpy as np
        p = pycoQC_plot(pycoQC_parse(summary_file=str(guppy_summary),
                                     quiet=True), sample=0, quiet=True)
        trace = p.channels_activity()["data"][0]
        assert np.asarray(trace["z"]).shape[0] == len(trace["y"])

    def test_channel_heatmap_counts_are_exact(self, guppy_summary):
        """The grid started from ones, so every empty cell reported a read
        that was never there and the totals were inflated by the grid size."""
        import numpy as np
        parser = pycoQC_parse(summary_file=str(guppy_summary), quiet=True)
        p = pycoQC_plot(parser, sample=0, quiet=True)
        helper = getattr(p, "_pycoQC_plot__channels_activity_data")

        _, reads = helper(df_level="all", count_level="reads", smooth_sigma=0)
        assert np.asarray(reads["z"][0]).sum() == len(parser.reads_df)

        _, bases = helper(df_level="all", count_level="bases", smooth_sigma=0)
        assert np.asarray(bases["z"][0]).sum() == parser.reads_df["read_len"].sum()


class TestPathTokenMatching:

    @pytest.mark.parametrize("path, token, expected", [
        ("/data/fastq_fail/r.fastq", "fail", True),
        ("/data/fastq_pass/r.fastq", "pass", True),
        ("/mnt/failsafe_backup/fastq_pass/r.fastq", "fail", False),
        ("/mnt/failover/fastq_pass/r.fastq", "fail", False),
        ("/data/passage_3/fastq_fail/r.fastq", "pass", False),
        ("FAW88668_pass_barcode06_0.pod5", "pass", True),
        ("FAW88668_fail_barcode06_0.pod5", "fail", True),
    ])
    def test_only_whole_words_match(self, path, token, expected):
        from pycoQC.common import path_has_token
        assert path_has_token(path, token) is expected

    def test_a_run_under_a_failover_mount_is_not_marked_failing(self, tmp_path):
        """Regression: the substring test marked every read of a good run as
        failing when the run happened to live under such a directory."""
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        d = tmp_path / "failover_mount" / "fastq_pass"
        d.mkdir(parents=True)
        with gzip.open(d / "r.fastq.gz", "wt") as fh:
            for i in range(5):
                fh.write(f"@read{i} runid=r1 ch=1 "
                         f"start_time=2026-01-01T00:0{i}:00Z\n"
                         "ACGTACGTAC\n+\nIIIIIIIIII\n")
        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(fastq_dir=str(d), seq_summary_fn=str(out),
                             threads=3)
        df = pd.read_csv(out, sep="\t")
        assert set(df["passes_filtering"]) == {True}


class TestFallbacksAreVisible:

    def test_unreadable_config_warns(self, tmp_path, guppy_summary, caplog):
        """A typo in --config_file silently produced a different report."""
        from pycoQC.pycoQC import pycoQC
        with caplog.at_level("WARNING"):
            pycoQC(summary_file=str(guppy_summary),
                   html_outfile=str(tmp_path / "r.html"),
                   config_file=str(tmp_path / "no_such_config.json"),
                   quiet=True)
        assert "configuration file" in caplog.text.lower()

    def test_unreadable_template_warns(self, tmp_path, guppy_summary, caplog):
        from pycoQC.pycoQC import pycoQC
        with caplog.at_level("WARNING"):
            pycoQC(summary_file=str(guppy_summary),
                   html_outfile=str(tmp_path / "r.html"),
                   template_file=str(tmp_path / "no_such_template.j2"),
                   quiet=True)
        assert "template file" in caplog.text.lower()


class TestCoverageSpansWholeAlignments:
    """An alignment used to be charged entirely to the bin its first base
    fell in, which turned real coverage on a small reference into spikes."""

    @staticmethod
    def _truth(starts, ends, total_len, nbins):
        import numpy as np
        depth = np.zeros(int(total_len))
        for s, e in zip(starts.astype(int), ends.astype(int)):
            depth[s:e] += 1
        edges = (np.arange(nbins + 1) * total_len / nbins).astype(int)
        return np.array([depth[edges[i]:edges[i + 1]].mean()
                         for i in range(nbins)])

    def test_matches_per_base_coverage(self):
        import numpy as np
        rng = np.random.default_rng(0)
        total_len, nbins = 30000, 500
        sites = np.linspace(0, total_len - 1000, 30).astype(int)
        starts = rng.choice(sites, 1000).astype(float)
        ends = starts + 1000
        got = pycoQC_plot._binned_coverage(starts, ends, total_len, nbins)
        np.testing.assert_allclose(got, self._truth(starts, ends, total_len,
                                                    nbins), atol=1e-9)

    def test_amplicon_reference_is_not_all_spikes(self):
        """The failure this replaces: 471 of 500 bins empty, peaks 17x the
        true depth, while the overall coverage line still looked right."""
        import numpy as np
        rng = np.random.default_rng(0)
        total_len, nbins = 30000, 500
        sites = np.linspace(0, total_len - 1000, 30).astype(int)
        starts = rng.choice(sites, 1000).astype(float)
        got = pycoQC_plot._binned_coverage(starts, starts + 1000, total_len,
                                           nbins)
        assert (got > 0).sum() > nbins * 0.9
        assert got.max() < 60

    def test_total_bases_are_conserved(self):
        import numpy as np
        rng = np.random.default_rng(1)
        total_len, nbins = 100000, 250
        starts = rng.integers(0, total_len - 500, 400).astype(float)
        ends = starts + rng.integers(100, 500, 400)
        got = pycoQC_plot._binned_coverage(starts, ends, total_len, nbins)
        assert got.sum() * (total_len / nbins) == pytest.approx(
            float((ends - starts).sum()))

    def test_alignment_inside_a_single_bin(self):
        import numpy as np
        got = pycoQC_plot._binned_coverage(np.array([0.0]), np.array([10.0]),
                                           1000, 10)
        assert got[0] == pytest.approx(0.1)
        assert got[1:].sum() == 0

    def test_alignment_spanning_everything(self):
        import numpy as np
        got = pycoQC_plot._binned_coverage(np.array([0.0]), np.array([1000.0]),
                                           1000, 10)
        np.testing.assert_allclose(got, 1.0)

    def test_short_reference_bins_track_the_labels(self):
        """Binning floor-divided a truncated step while the labels used the
        true width, so on a short reference the two drifted apart and the
        last bin absorbed everything past the end."""
        import numpy as np
        total_len, nbins = 1200, 500
        starts = np.array([600.0])
        got = pycoQC_plot._binned_coverage(starts, starts + 1, total_len, nbins)
        # The label for position 600 sits at 600 * nbins / total_len = 250.
        assert int(np.argmax(got)) == 250


class TestFast5FileTypeAttribute:

    def test_multi_read_written_as_str(self, tmp_path):
        """ont_fast5_api writes the attribute as str, which failed a
        comparison against bytes and sent the file down the single-read path
        until it died looking for /Raw/Reads."""
        from conftest import _write_basecalled_fast5
        from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
        h5py = pytest.importorskip("h5py")

        d = tmp_path / "f5"
        d.mkdir()
        path = _write_basecalled_fast5(d / "batch.fast5", n_reads=5)
        with h5py.File(path, "r+") as fh:
            fh.attrs["file_type"] = "multi-read"      # a plain str

        out = tmp_path / "s.txt"
        Fast5_to_seq_summary(fast5_dir=str(d), seq_summary_fn=str(out),
                             threads=3)
        assert len(pd.read_csv(out, sep="\t")) == 5


class TestBarcodePlaceholders:

    def test_bam_without_barcodes_has_no_barcode(self, tmp_path):
        """Writing "none" made an unbarcoded run look like a barcoded one
        with a single sample called none."""
        from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
        pysam = pytest.importorskip("pysam")
        d = tmp_path / "bam"
        d.mkdir()
        header = {"HD": {"VN": "1.6", "SO": "unknown"},
                  "RG": [{"ID": "rg1", "DS": "runid=r1 basecall_model=m1"}]}
        with pysam.AlignmentFile(str(d / "r.bam"), "wb", header=header) as bam:
            for i in range(3):
                rec = pysam.AlignedSegment()
                rec.query_name = f"read{i}"
                rec.query_sequence = "ACGT" * 25
                rec.query_qualities = pysam.qualitystring_to_array("I" * 100)
                rec.flag = 4
                rec.set_tags([("RG", "rg1", "Z"), ("ch", 1, "i"),
                              ("st", "2026-01-01T00:00:00Z", "Z")])
                bam.write(rec)
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert df["barcode_arrangement"].isna().all()

    def test_duplicate_ids_in_a_barcode_file_do_not_multiply_reads(self, tmp_path):
        """A left merge on a duplicated key repeats the read instead of
        failing, so every count came out high."""
        summary = write_summary(tmp_path / "s.txt", n=600, barcodes=False)
        barcodes = tmp_path / "bc.txt"
        ids = [f"r{i:06d}" for i in range(600)] + ["r000001", "r000002"]
        pd.DataFrame({
            "read_id": ids,
            "barcode_arrangement": ["barcode01"] * len(ids),
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        assert len(parser.reads_df) == 600


class TestBareOutputFilename:

    def test_converter_accepts_a_name_with_no_directory(self, tmp_path,
                                                        fastq_dir, monkeypatch):
        """dirname of a bare name is "", which is never writable."""
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        monkeypatch.chdir(tmp_path)
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir),
                             seq_summary_fn="summary.txt", threads=3)
        assert (tmp_path / "summary.txt").exists()


class TestSmoothingKeepsMeasurements:

    @pytest.mark.parametrize("sigma", [1.0, 1.5, 2.0, 3.0])
    def test_an_isolated_measured_bin_survives(self, sigma):
        """The weight threshold governs how far a value carries into a gap;
        on its own it discarded lone measurements at higher sigma."""
        import numpy as np
        values = np.full(200, np.nan)
        values[:60] = 5.0
        values[150] = 9.0           # isolated, far from the dense stretch
        smoothed = pycoQC_plot._smooth_with_gaps(values, sigma=sigma)
        assert np.isfinite(smoothed[150])

# -*- coding: utf-8 -*-
"""Parsing of what the basecallers actually write.

Guppy, Dorado and MinKNOW disagree about fastq headers and BAM tags, and the
Oxford Nanopore output specification differs from what Dorado emits. These
pin the shapes seen in the wild so a parser change cannot quietly drop one.
"""

import gzip
import random

import pandas as pd
import pytest

from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
from pycoQC.Fastq_to_seq_summary import FastqParser, Fastq_to_seq_summary


def df_bases(summary):
    return int(pd.read_csv(summary, sep="\t")["sequence_length_template"].sum())


def one_record(header, seq="ACGTACGTAC", qual="IIIIIIIIII"):
    return FastqParser.read_entry(
        [header + "\n", seq + "\n", "+\n", qual + "\n"], "pass")


class TestFastqHeaderDialects:

    def test_guppy_key_value(self):
        parsed = one_record(
            "@read1 runid=acfd15dc sampleid=s1 read=2743 ch=207 "
            "start_time=2022-06-27T23:14:28Z barcode=barcode17")
        d = parsed["read1"]
        assert d["run_id"] == "acfd15dc"
        assert d["channel"] == "207"
        assert d["barcode_arrangement"] == "barcode17"
        assert d["start_time"].year == 2022

    def test_dorado_sam_style_tags(self):
        """Regression: splitting on "=" raised ValueError on the first read,
        so Dorado and current MinKNOW output could not be parsed at all."""
        parsed = one_record(
            "@bd8655fb RG:Z:e4994c62_dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "
            "ch:i:123 st:Z:2025-01-06T10:06:36.778368+00:00 du:f:1.23 "
            "qs:f:12.5 PU:Z:FAW88668 SM:Z:barcode02")
        d = parsed["bd8655fb"]
        assert d["channel"] == "123"
        # Dorado has no runid field; it is the first part of the read group.
        assert d["run_id"] == "e4994c62"
        assert d["barcode_arrangement"] == "barcode02"
        assert d["start_time"].year == 2025

    def test_dorado_tab_separated(self):
        parsed = one_record(
            "@bd8655fb\tch:i:123\tdu:f:1.23\t"
            "st:Z:2025-01-06T10:06:36.778368+00:00\tal:Z:barcode07")
        d = parsed["bd8655fb"]
        assert d["channel"] == "123"
        assert d["barcode_arrangement"] == "barcode07"

    def test_illumina_casava_header(self):
        """Illumina carries none of the nanopore fields, and its
        "1:N:0:INDEX" block must not be mistaken for a tag."""
        parsed = one_record(
            "@A00111:222:HXXXXDSXX:1:1101:1000:2000 1:N:0:ATCGATCG")
        assert "A00111:222:HXXXXDSXX:1:1101:1000:2000" in parsed
        d = parsed["A00111:222:HXXXXDSXX:1:1101:1000:2000"]
        assert d["run_id"] == "" and d["channel"] == "" and d["start_time"] == ""
        assert d["sequence_length_template"] == 10

    def test_plain_header_without_metadata(self):
        parsed = one_record("@plainread")
        assert parsed["plainread"]["sequence_length_template"] == 10

    def test_unknown_fields_are_ignored(self):
        parsed = one_record("@read4 runid=abc ch=5 odd=x:y:z bare:token")
        d = parsed["read4"]
        assert d["run_id"] == "abc" and d["channel"] == "5"

    @pytest.mark.parametrize("barcode_field, expected", [
        ("barcode=barcode01", "barcode01"),
        ("SM:Z:barcode02", "barcode02"),
        ("BC:Z:barcode03", "barcode03"),
        ("al:Z:barcode04", "barcode04"),
    ])
    def test_barcode_read_from_any_dialect(self, barcode_field, expected):
        parsed = one_record("@r ch=1 " + barcode_field)
        assert parsed["r"]["barcode_arrangement"] == expected


class TestParseHeaderFields:

    def test_value_containing_equals_is_kept_whole(self):
        fields = FastqParser.parse_header_fields(
            ["RG:Z:run_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"])
        assert fields["RG"] == "run_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

    def test_value_containing_colons_is_kept_whole(self):
        fields = FastqParser.parse_header_fields(
            ["st:Z:2025-01-06T10:06:36.778368+00:00"])
        assert fields["st"] == "2025-01-06T10:06:36.778368+00:00"

    @pytest.mark.parametrize("item", ["", "nocolonnoequals", "a:b", ":::"])
    def test_malformed_items_do_not_raise(self, item):
        FastqParser.parse_header_fields([item])


class TestBamReadGroups:

    def _header(self, **overrides):
        rg = {
            "ID": "e2419f12_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_barcode03",
            "DS": "runid=e2419f12 basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
            "PU": "FAW88668", "PM": "promethion", "LB": "sample1",
            "SM": "barcode03", "al": "alias1",
        }
        rg.update(overrides)
        pysam = pytest.importorskip("pysam")
        return pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6"}, "RG": [rg]})

    def test_model_is_the_model_not_the_run_id(self):
        """Regression: DS was read positionally and Oxford Nanopore list
        runid first, so the model column received the run id."""
        groups = Bam_to_seq_summary.parse_read_groups(self._header())
        group = next(iter(groups.values()))
        assert group["dorado_model"] == "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
        assert group["run_id"] == "e2419f12"

    def test_barcode_from_read_group(self):
        """MinKNOW writes no BC tag; the barcode is SM on the read group."""
        groups = Bam_to_seq_summary.parse_read_groups(self._header())
        assert next(iter(groups.values()))["barcode"] == "barcode03"

    def test_alias_used_when_no_sm(self):
        header = self._header()
        groups = Bam_to_seq_summary.parse_read_groups(header)
        assert next(iter(groups.values()))["barcode"] == "barcode03"

    def test_missing_fields_give_blanks_not_errors(self):
        """A BAM from another tool may carry almost nothing."""
        pysam = pytest.importorskip("pysam")
        header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6"}, "RG": [{"ID": "bare"}]})
        group = Bam_to_seq_summary.parse_read_groups(header)["bare"]
        assert group["flowcell_id"] == "" and group["device_id"] == ""
        assert group["dorado_model"] == "" and group["barcode"] == ""
        assert group["run_id"] == "bare"

    def test_no_read_groups_at_all(self):
        pysam = pytest.importorskip("pysam")
        header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
        assert Bam_to_seq_summary.parse_read_groups(header) == {}


@pytest.fixture(scope="module")
def illumina_summary(tmp_path_factory):
    """A summary built from Illumina reads, which carry no run metadata."""
    d = tmp_path_factory.mktemp("illumina")
    rng = random.Random(11)
    lengths = []
    with gzip.open(d / "sample_R1.fastq.gz", "wt") as fh:
        for i in range(600):
            # Trimmed Illumina reads vary in length below the cycle count.
            n = rng.randint(90, 151)
            lengths.append(n)
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            qual = "".join(chr(rng.randint(25, 40) + 33) for _ in range(n))
            fh.write(f"@A00111:222:HXXXXDSXX:1:1101:{1000+i}:{2000+i}"
                     f" 1:N:0:ATCGATCG\n{seq}\n+\n{qual}\n")
    out = d / "summary.txt"
    Fastq_to_seq_summary(fastq_dir=str(d), seq_summary_fn=str(out), threads=3)
    return out


class TestUniversalFastq:

    def test_summary_has_empty_run_metadata(self, illumina_summary):
        df = pd.read_csv(illumina_summary, sep="\t")
        assert len(df) == 600
        for column in ("run_id", "channel", "start_time"):
            assert df[column].isna().all()
        assert df["sequence_length_template"].between(90, 151).all()

    def test_report_builds_without_run_metadata(self, illumina_summary, tmp_path):
        """Regression: every read was discarded by the NA filter, because
        run_id, channel and start_time were required."""
        from pycoQC.pycoQC import pycoQC
        import json

        html = tmp_path / "report.html"
        js = tmp_path / "report.json"
        pycoQC(summary_file=str(illumina_summary), html_outfile=str(html),
               json_outfile=str(js), quiet=True)

        assert html.exists() and html.stat().st_size > 50_000
        basecall = json.loads(js.read_text())["All Reads"]["basecall"]
        assert basecall["reads_number"] == 600
        assert 90 <= basecall["N50"] <= 151
        assert basecall["bases_number"] == df_bases(illumina_summary)

    def test_time_and_channel_plots_are_skipped(self, illumina_summary):
        from pycoQC.common import pycoQCError
        from pycoQC.pycoQC_parse import pycoQC_parse
        from pycoQC.pycoQC_plot import pycoQC_plot

        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(illumina_summary),
                                           quiet=True), quiet=True)
        for method in ("output_over_time", "read_len_over_time",
                       "channels_activity"):
            with pytest.raises(pycoQCError):
                getattr(plotter, method)()

    def test_length_and_quality_plots_still_work(self, illumina_summary):
        from pycoQC.pycoQC_parse import pycoQC_parse
        from pycoQC.pycoQC_plot import pycoQC_plot

        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(illumina_summary),
                                           quiet=True), quiet=True)
        for method in ("read_len_1D", "read_qual_1D", "basecall_summary"):
            assert getattr(plotter, method)() is not None

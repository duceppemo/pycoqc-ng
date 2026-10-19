# -*- coding: utf-8 -*-
"""Splitting a sequencing summary into one file per barcode."""

import pandas as pd
import pytest

from pycoQC.Barcode_split import barcode_split
from pycoQC.common import pycoQCError


def split_files(output_dir):
    return sorted(p.name for p in output_dir.iterdir() if p.is_file())


class TestBarcodeSplit:

    def test_writes_one_file_per_barcode(self, tmp_path, barcoded_summary):
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), quiet=True)
        names = split_files(out)
        assert "sequencing_summary_barcode01.txt" in names
        assert "sequencing_summary_barcode02.txt" in names

    def test_rows_are_partitioned_without_loss(self, tmp_path, barcoded_summary):
        """Every classified read must land in exactly one output file."""
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), output_unclassified=True,
                      min_barcode_percent=0, quiet=True)

        source = pd.read_csv(barcoded_summary, sep="\t")
        written = pd.concat([pd.read_csv(p, sep="\t")
                             for p in out.iterdir() if p.is_file()])
        assert len(written) == len(source)
        assert set(written["read_id"]) == set(source["read_id"])

    def test_each_file_holds_exactly_its_own_reads(self, tmp_path, barcoded_summary):
        """Barcode columns are dropped on the way out, since the filename
        already carries the barcode, so check the membership instead."""
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), output_unclassified=True,
                      min_barcode_percent=0, quiet=True)

        source = pd.read_csv(barcoded_summary, sep="\t")
        for path in sorted(p for p in out.iterdir() if p.is_file()):
            barcode = path.stem.replace("sequencing_summary_", "")
            df = pd.read_csv(path, sep="\t")
            assert not [c for c in df.columns if c.startswith("barcode")]
            expected = source.loc[source["barcode_arrangement"] == barcode,
                                  "read_id"]
            assert set(df["read_id"]) == set(expected)

    def test_unclassified_skipped_by_default(self, tmp_path, barcoded_summary):
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), min_barcode_percent=0, quiet=True)
        assert not any("unclassified" in n for n in split_files(out))

    def test_unclassified_written_when_asked(self, tmp_path, barcoded_summary):
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), output_unclassified=True,
                      min_barcode_percent=0, quiet=True)
        assert any("unclassified" in n for n in split_files(out))

    def test_low_frequency_barcodes_are_dropped(self, tmp_path, barcoded_summary):
        """barcode03 holds ~11% of reads, so a 20% cutoff must exclude it
        while the two larger barcodes survive."""
        out = tmp_path / "split"
        barcode_split(summary_file=str(barcoded_summary),
                      output_dir=str(out), min_barcode_percent=20, quiet=True)
        names = split_files(out)
        assert "sequencing_summary_barcode01.txt" in names
        assert "sequencing_summary_barcode03.txt" not in names

    def test_summary_without_barcodes_is_rejected(self, tmp_path, guppy_summary):
        with pytest.raises(pycoQCError, match="[Nn]o barcode"):
            barcode_split(summary_file=str(guppy_summary),
                          output_dir=str(tmp_path / "split"), quiet=True)

# -*- coding: utf-8 -*-
"""The shared helpers: argument checking, file discovery and merging."""

import gzip
import logging

import pandas as pd
import pytest

from pycoQC.common import (check_arg, dict_to_str, expand_file_names,
                           get_logger, head, is_readable_file,
                           merge_files_to_df, mkbasedir, pycoQCError,
                           recursive_file_gen, sequencing_summary_file_sample)


def write_summary(path, n=50, run_id="run1", gzipped=False):
    df = pd.DataFrame({
        "read_id": [f"{run_id}_r{i}" for i in range(n)],
        "run_id": run_id,
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
        "sequence_length_template": [100 + i for i in range(n)],
        "mean_qscore_template": [10.0 + (i % 5) for i in range(n)],
    })
    opener = gzip.open if gzipped else open
    with opener(path, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)
    return path


class TestCheckArg:

    def test_returns_a_value_of_the_right_type(self):
        assert check_arg("n", 5, int) == 5
        assert check_arg("s", "x", str) == "x"

    def test_coerces_when_it_can(self):
        assert check_arg("n", "5", int) == 5
        assert check_arg("f", 5, float) == 5.0

    def test_rejects_what_it_cannot_coerce(self):
        with pytest.raises(Exception, match="cannot be coerced"):
            check_arg("n", "not a number", int)

    def test_none_allowed_by_default(self):
        assert check_arg("n", None, int) is None

    def test_none_rejected_when_not_allowed(self):
        with pytest.raises(Exception):
            check_arg("n", None, int, allow_none=False)

    def test_minimum_is_enforced(self):
        assert check_arg("n", 5, int, min=1) == 5
        with pytest.raises(Exception, match="too low"):
            check_arg("n", 0.5, float, min=1)

    def test_maximum_is_enforced(self):
        with pytest.raises(Exception, match="too high"):
            check_arg("n", 100, int, max=60)

    def test_choices_are_enforced(self):
        assert check_arg("c", "a", str, choices=["a", "b"]) == "a"
        with pytest.raises(Exception, match="choices"):
            check_arg("c", "z", str, choices=["a", "b"])

    def test_zero_passes_a_minimum_of_zero(self):
        """A min of 0 is falsy, so the bound is not applied; 0 is valid."""
        assert check_arg("n", 0, int, min=0) == 0


class TestFileDiscovery:

    def test_is_readable_file(self, tmp_path):
        existing = write_summary(tmp_path / "a.txt")
        assert is_readable_file(str(existing))
        assert not is_readable_file(str(tmp_path / "missing.txt"))
        assert not is_readable_file(str(tmp_path))

    def test_expand_a_single_name(self, tmp_path):
        path = write_summary(tmp_path / "a.txt")
        assert expand_file_names(str(path)) == [str(path)]

    def test_expand_a_glob(self, tmp_path):
        write_summary(tmp_path / "a.txt")
        write_summary(tmp_path / "b.txt")
        found = expand_file_names(str(tmp_path / "*.txt"))
        assert len(found) == 2

    def test_expand_a_list_of_names(self, tmp_path):
        a = write_summary(tmp_path / "a.txt")
        b = write_summary(tmp_path / "b.txt")
        assert sorted(expand_file_names([str(a), str(b)])) == sorted(
            [str(a), str(b)])

    def test_expand_a_list_holding_one_glob(self, tmp_path):
        write_summary(tmp_path / "a.txt")
        write_summary(tmp_path / "b.txt")
        assert len(expand_file_names([str(tmp_path / "*.txt")])) == 2

    def test_missing_file_is_reported(self, tmp_path):
        with pytest.raises(pycoQCError, match="No files found"):
            expand_file_names(str(tmp_path / "nothing_here.txt"))

    def test_wrong_type_is_reported(self):
        with pytest.raises(pycoQCError, match="has to be"):
            expand_file_names(42)


class TestRecursiveFileGen:

    def test_finds_files_at_the_top_level(self, tmp_path):
        (tmp_path / "a.fast5").touch()
        (tmp_path / "b.fast5").touch()
        assert len(list(recursive_file_gen(str(tmp_path), "fast5"))) == 2

    def test_descends_until_it_finds_matches(self, tmp_path):
        deep = tmp_path / "x" / "y" / "z"
        deep.mkdir(parents=True)
        (deep / "a.fast5").touch()
        assert len(list(recursive_file_gen(str(tmp_path), "fast5"))) == 1

    def test_accepts_several_extensions(self, tmp_path):
        (tmp_path / "a.fastq").touch()
        (tmp_path / "b.fq.gz").touch()
        found = list(recursive_file_gen(
            str(tmp_path), ["fq", "fq.gz", "fastq", "fastq.gz"]))
        assert len(found) == 2

    def test_other_extensions_are_ignored(self, tmp_path):
        (tmp_path / "a.txt").touch()
        assert list(recursive_file_gen(str(tmp_path), "fast5")) == []


class TestMergeFilesToDf:

    def test_single_file(self, tmp_path):
        path = write_summary(tmp_path / "a.txt", n=10)
        assert len(merge_files_to_df([str(path)])) == 10

    def test_several_files_are_concatenated(self, tmp_path):
        a = write_summary(tmp_path / "a.txt", n=10, run_id="r1")
        b = write_summary(tmp_path / "b.txt", n=15, run_id="r2")
        df = merge_files_to_df([str(a), str(b)])
        assert len(df) == 25
        assert set(df["run_id"]) == {"r1", "r2"}

    def test_merge_keeps_only_shared_columns(self, tmp_path):
        a = write_summary(tmp_path / "a.txt", n=5)
        extra = pd.read_csv(a, sep="\t")
        extra["only_here"] = 1
        b = tmp_path / "b.txt"
        extra.to_csv(b, sep="\t", index=False)
        df = merge_files_to_df([str(a), str(b)])
        assert "only_here" not in df.columns

    def test_empty_file_is_reported(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("read_id\trun_id\n")
        with pytest.raises(pycoQCError, match="No valid read"):
            merge_files_to_df([str(empty)])


class TestSummarySampling:

    def test_returns_a_dataframe_when_no_outfile(self, tmp_path, capsys):
        path = write_summary(tmp_path / "a.txt", n=100)
        df = sequencing_summary_file_sample(str(path), n_seq=10)
        assert len(df) == 10

    def test_writes_a_file_when_asked(self, tmp_path, capsys):
        path = write_summary(tmp_path / "a.txt", n=100)
        out = tmp_path / "sampled.txt"
        sequencing_summary_file_sample(str(path), outfile=str(out), n_seq=20)
        assert len(pd.read_csv(out, sep="\t")) == 20

    def test_writes_gzip_when_the_name_says_so(self, tmp_path, capsys):
        path = write_summary(tmp_path / "a.txt", n=100)
        out = tmp_path / "sampled.txt.gz"
        sequencing_summary_file_sample(str(path), outfile=str(out), n_seq=20)
        with gzip.open(out, "rt") as fh:
            assert len(pd.read_csv(fh, sep="\t")) == 20

    def test_samples_each_run_proportionally(self, tmp_path, capsys):
        a = write_summary(tmp_path / "a.txt", n=90, run_id="r1")
        extra = write_summary(tmp_path / "b.txt", n=10, run_id="r2")
        both = pd.concat([pd.read_csv(a, sep="\t"), pd.read_csv(extra, sep="\t")])
        merged = tmp_path / "both.txt"
        both.to_csv(merged, sep="\t", index=False)

        df = sequencing_summary_file_sample(str(merged), n_seq=50)
        counts = df["run_id"].value_counts()
        assert counts["r1"] > counts["r2"]


class TestSmallHelpers:

    def test_dict_to_str_formats_pairs(self):
        out = dict_to_str({"a": 1, "b": 2})
        assert "a" in out and "1" in out and "b" in out

    def test_dict_to_str_handles_empty(self):
        assert isinstance(dict_to_str({}), str)

    def test_head_returns_the_first_lines(self, tmp_path, capsys):
        path = write_summary(tmp_path / "a.txt", n=50)
        head(str(path), n=5)
        assert capsys.readouterr().out.count("\n") >= 5

    def test_mkbasedir_creates_the_parent(self, tmp_path):
        target = tmp_path / "deep" / "deeper" / "file.txt"
        mkbasedir(str(target), exist_ok=True)
        assert target.parent.is_dir()

    def test_mkbasedir_tolerates_a_bare_name(self):
        """No directory part means there is nothing to create."""
        mkbasedir("bare_name.txt", exist_ok=True)

    @pytest.mark.parametrize("verbose, quiet, expected", [
        (True, False, logging.DEBUG),
        (False, True, logging.WARNING),
        (False, False, logging.INFO),
    ])
    def test_logger_level_follows_the_flags(self, verbose, quiet, expected):
        logger = get_logger(name=f"t{verbose}{quiet}", verbose=verbose,
                            quiet=quiet)
        assert logger.level == expected


class TestDocumentationHelpers:
    """doc_func, make_arg_dict and arg_opt drive both the notebook help and
    the construction of the command line parsers."""

    @staticmethod
    def sample(first: int = 5, flag: bool = False, names: list = (),
               text: str = "x"):
        """One line of description.

        A second line, still description.

        * first
            The first argument.
        * flag
            A boolean argument.
        * names
            A list argument.
        * text
            A string argument.
        """

    def test_doc_func_stops_at_the_argument_list(self):
        from pycoQC.common import doc_func
        described = doc_func(self.sample)
        assert "One line of description." in described
        assert "still description" in described
        assert "The first argument" not in described

    def test_make_arg_dict_reads_defaults_and_types(self):
        from pycoQC.common import make_arg_dict
        args = make_arg_dict(self.sample)
        assert args["first"]["default"] == 5
        assert args["first"]["type"] is int
        assert "The first argument" in args["first"]["help"]

    def test_make_arg_dict_accepts_a_class(self):
        from pycoQC.common import make_arg_dict
        from pycoQC.pycoQC_parse import pycoQC_parse
        args = make_arg_dict(pycoQC_parse)
        assert "summary_file" in args

    def test_arg_opt_turns_a_false_flag_into_store_true(self):
        from pycoQC.common import arg_opt
        opt = arg_opt(self.sample, "flag")
        assert opt["action"] == "store_true"
        assert "type" not in opt

    def test_arg_opt_turns_a_list_into_nargs(self):
        from pycoQC.common import arg_opt
        opt = arg_opt(self.sample, "names")
        assert opt["nargs"] == "*"

    def test_arg_opt_appends_the_default_to_the_help(self):
        from pycoQC.common import arg_opt
        opt = arg_opt(self.sample, "first")
        assert "%(default)s" in opt["help"]

    def test_jhelp_renders_without_a_display(self, monkeypatch):
        """jhelp writes Markdown through IPython, which is an optional extra;
        it should not be imported unless jhelp is called."""
        pytest.importorskip("IPython")
        from pycoQC.common import jhelp
        shown = []
        import IPython.display
        monkeypatch.setattr(IPython.display, "display",
                            lambda *a, **k: shown.append(a))
        jhelp(self.sample)
        assert shown

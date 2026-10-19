# Development

## Setup

```bash
git clone https://github.com/duceppemo/pycoqc-ng.git
cd pycoqc-ng
pip install -e ".[dev,notebook]"
```

Or with conda, which saves compiling `pysam`:

```bash
conda create -n pycoqc-ng-dev -c conda-forge -c bioconda \
  python=3.11 pysam pod5 pandas numpy scipy plotly jinja2 h5py tqdm python-dateutil pip
conda activate pycoqc-ng-dev
pip install -e ".[dev,notebook]"
```

## Tests

```bash
pytest                      # whole suite
pytest tests/test_cli.py    # one file
pytest -k barcode           # by name
pytest --cov --cov-report=term
```

Fixtures in `tests/conftest.py` synthesise FASTQ and Dorado-style BAM so the suite needs
no large inputs. Tests that need a real summary use the small one under
`docs/pycoQC/data/`, and skip if it is absent.

Coverage is configured for multiprocessing in `pyproject.toml`. The converters parse in
worker processes, so without `concurrency = ["multiprocessing", "thread"]` that code
reads as untested even when it runs.

## Linting

```bash
python -m pyflakes pycoQC tests
```

This runs in CI and must stay clean. It only became useful once the star imports were
removed: pyflakes cannot resolve names coming from `import *`, so prefer explicit
imports.

## Layout

```
pycoQC/
  __main__.py              argument parsing for all six commands
  pycoQC.py                the top-level entry point
  pycoQC_parse.py          reads summary, barcode and BAM files into dataframes
  pycoQC_plot.py           one method per plot, each returning a plotly figure
  pycoQC_report.py         assembles plots into HTML and JSON
  common.py                shared helpers, errors and logging
  *_to_seq_summary.py      the four converters
  Barcode_split.py         splits a summary per barcode
  templates/               default config and Jinja2 template
recipe/                    conda recipe, mirrored into bioconda
docs/wiki/                 these pages
tests/
```

Adding a plot means a method on `pycoQC_plot` returning a plotly figure, plus an entry
in `templates/pycoQC_config.json`.

## CI

* `ci.yml` runs the suite on Python 3.9 to 3.13 plus macOS, lints with pyflakes, builds
  the distributions and checks that a clean install from the sdist works outside the
  repository.
* `publish.yml` publishes to PyPI via Trusted Publishing. A manual dispatch can target
  TestPyPI first.
* `wiki.yml` syncs `docs/wiki/` to the GitHub wiki on push.

## Releasing

1. Update `__version__` in `pycoQC/__init__.py`, which `pyproject.toml` reads.
2. Add a `CHANGELOG.md` entry.
3. Update `version`, `date-released` and the DOI list in `CITATION.cff`.
4. Tag and push: `git tag -a vX.Y.Z -m "..." && git push origin vX.Y.Z`.
5. Publish a GitHub release for the tag. That uploads to PyPI and, through the Zenodo
   webhook, archives the release and mints a DOI.
6. Update `recipe/meta.yaml`: new `version`, `build: number: 0`, and the sha256 of the
   new sdist:

   ```bash
   curl -s https://pypi.org/pypi/pycoqc-ng/json \
     | python -c "import sys,json; print([f['digests']['sha256'] for f in json.load(sys.stdin)['urls'] if f['packagetype']=='sdist'][0])"
   ```

7. Open a PR against [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
   updating `recipes/pycoqc-ng/meta.yaml`.

Zenodo only archives releases published *after* its webhook was enabled. If a release
is missed, deleting and recreating it on GitHub re-fires the event; the publish workflow
passes `skip-existing`, so the repeated PyPI upload is a no-op.

## Conventions

Match the surrounding code rather than reformatting it. Explicit imports, no `import *`.
Comments should explain why something is done, especially where it guards against a bug
that has already happened once.

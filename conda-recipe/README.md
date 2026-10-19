# Conda recipe

This is the recipe for the bioconda package, kept here so it can be built and
tested locally:

```sh
conda build conda-recipe -c conda-forge -c bioconda
```

Bioconda builds from its own copy at `recipes/pycoqc-ng/meta.yaml` in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes), which is the
source of truth once the package is published. Keep the two in step when
releasing: bump `version`, reset `build: number` to 0, and replace `sha256`
with the hash of the new sdist on PyPI:

```sh
curl -s https://pypi.org/pypi/pycoqc-ng/json \
  | python -c "import sys,json; print(json.load(sys.stdin)['urls'][1]['digests']['sha256'])"
```

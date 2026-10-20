# Configuration

Which plots the HTML report contains, and how each is drawn, comes from a JSON
configuration file. Without one, a built-in default is used: one for a basecalled run,
and a shorter one for a run that has not been basecalled, where most of the usual plots
have nothing to draw.

## Starting from the default

```bash
pycoQC --default_config > my_config.json
pycoQC -f sequencing_summary.txt -o report.html --config_file my_config.json
```

If `--config_file` is not given, `pycoQC` looks for `~/.pycoQC` and then
`~/.config/pycoQC/config`, falling back to the built-in default. Dropping a file at one
of those paths changes the default for every run on that machine.

**A configuration you supply is used as it stands**, whether you pass it with
`--config_file` or leave it in your home directory. It therefore also replaces the
shorter default that a run without basecalls would otherwise get, and such a run will
then show a report of mostly skipped plots. If you work with POD5 or raw FAST5, either
leave the default in place for those runs or keep a second configuration for them; the
built-in one is `pycoQC_config_signal.json` inside the installed package, listing
`run_summary`, `signal_summary`, `read_duration_1D`, `output_over_time`,
`read_duration_over_time`, `channels_activity` and `barcode_counts`.

## Structure

First-level keys name the plots, in the order they appear in the report. Second-level
keys are arguments passed to the corresponding plotting function.

```json
{
  "run_summary": {
    "plot_title": "General run summary"
  },
  "read_len_1D": {
    "plot_title": "Basecalled reads length",
    "color": "lightsteelblue",
    "nbins": 200,
    "smooth_sigma": 1
  },
  "channels_activity": {
    "plot_title": "Channel activity over time",
    "colorscale": [[0.0,"rgba(255,255,255,0)"], [1.0,"rgb(30,30,70)"]],
    "smooth_sigma": 1
  }
}
```

Remove a key to drop that plot; reorder keys to reorder the report. See
[Outputs](Outputs) for the full catalogue of plot names.

Common parameters:

| Parameter | Meaning |
|---|---|
| `plot_title` | Heading shown above the plot |
| `color` / `colorscale` | Any CSS colour or plotly colour scale |
| `nbins` | Histogram bins |
| `smooth_sigma` | Gaussian smoothing; `0` disables it |
| `time_bins` | Buckets for the over-time plots |

A plot whose data is unavailable is skipped with a note in the log, so leaving
alignment plots configured without passing `-a` is harmless.

## Examples

**A short report, length and quality only**

```json
{
  "run_summary":   { "plot_title": "Run summary" },
  "basecall_summary": { "plot_title": "Basecall summary" },
  "read_len_1D":   { "plot_title": "Read length" },
  "read_qual_1D":  { "plot_title": "Read quality" }
}
```

**Less smoothing and finer bins**

```json
{
  "read_len_1D": { "plot_title": "Read length", "nbins": 500, "smooth_sigma": 0 }
}
```

To see which parameters a plot accepts, ask the plotting function itself:

```python
from pycoQC.common import jhelp
from pycoQC.pycoQC_plot import pycoQC_plot

jhelp(pycoQC_plot.read_len_1D)
```

## Templates

`--template_file` replaces the Jinja2 template used for the page. The default is
`pycoQC/templates/spectre.html.j2` in the installed package; copy it and adjust rather
than starting from scratch. This changes the page around the plots, not the plots
themselves.

## Performance

`--sample` controls how many reads the plots are built from, 100,000 by default and
sampled deterministically, so a report is reproducible. Summary statistics always use
every read, so sampling changes the plots' resolution, not the numbers. On a very large
run, lowering it is the quickest way to cut both runtime and file size; `--sample 0`
disables sampling entirely.

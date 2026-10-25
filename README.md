# AgriStatsPy

[![CI](https://github.com/agristatspy/agristatspy/actions/workflows/ci.yml/badge.svg)](https://github.com/agristatspy/agristatspy/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%20|%203.11%20|%203.12-blue)
![License](https://img.shields.io/badge/license-MIT-green)

**Experimental design and statistical analysis for agricultural field trials in Python.**

AgriStatsPy provides a consistent, well-documented API for generating randomized
field layouts and analyzing the resulting data. It is built for agronomy
students, researchers, and extension professionals who need reproducible,
statistically sound results without writing boilerplate for every experiment.

## Features

- **Design generation** — CRD, RCBD, Latin square, factorial (CRD or RCBD layout) and split-plot designs with reproducible, seeded randomization.
- **Field layouts** — printable text layouts with block separators, Matplotlib layout figures (PNG/SVG), and CSV/Excel export.
- **Descriptive statistics** — per-treatment summaries and coefficient of variation.
- **One-command analysis** — `ag.analyze` chains normality check → ANOVA → post-hoc letters → mean ± SEM, SEm±, CD and CV%, with thesis-ready export (Excel, CSV, Markdown, HTML, LaTeX).
- **Analysis of variance** — CRD, RCBD, Latin square, factorial and split-plot ANOVA (with correct split-plot error terms), CV%, and explicit model assumptions; `ag.anova(...)` dispatches by design.
- **Mean comparison** — Tukey HSD, Fisher's LSD, Duncan's DMRT and Bonferroni, using the design's own error term, with compact letter displays.
- **Diagnostics and plots** — residual normality/variance checks, mean, box, interaction and residual plots.
- **Reports and Excel export** — Markdown, HTML and multi-sheet Excel workbooks of designs, ANOVA, comparisons and diagnostics.
- **Strict validation** — invalid inputs and incomplete data raise clear, package-specific exceptions.

## Installation

Requires Python 3.10 or later.

```bash
pip install agristatspy

# Optional: Excel export
pip install "agristatspy[excel]"
```

To install from source:

```bash
git clone https://github.com/agristatspy/agristatspy.git
cd agristatspy
pip install -e ".[dev]"
```

See [docs/installation.md](docs/installation.md) for details.

## Quick start

### Generate a design

```python
import agristatspy as ag

# Completely Randomized Design
crd = ag.design.crd(treatments=["Control", "N50", "N100"], replications=4, seed=12)
print(crd.text_layout())

# Randomized Complete Block Design (replications = number of blocks)
rcbd = ag.design.rcbd(
    treatments=["Control", "N50", "N100", "N150"],
    replications=3,
    seed=42,
)
print(rcbd.text_layout())

rcbd.to_csv("field_layout.csv")
```

### Analyze data in one command

```python
import agristatspy as ag
from agristatspy.datasets import load_rcbd_example

df = load_rcbd_example()
report = ag.analyze(df, "yield", "rcbd", treatment="treatment", block="block")
print(report)                      # normality -> ANOVA -> post-hoc letters, CV%, SEm±, CD

report.thesis_table()              # Mean ± SEM with letters, ready for a paper
report.export("results.xlsx")      # also .csv, .md, .html (Word), .tex
```

Choose the mean-comparison test with `posthoc=`: `"tukey"` (default), `"lsd"`, `"dmrt"`, `"bonferroni"`, a list like `["lsd", "dmrt"]`, `"all"` (every test, letters side by side), `"ask"` (prompts you) or `None`.

```python
report = ag.analyze(df, "yield", "rcbd", treatment="treatment", block="block", posthoc="all")
report.comparisons["dmrt"]         # each test's full pairwise table
```

`analyze` runs Shapiro–Wilk on the residuals first. If they are not normal it asks whether to proceed or review the data (in scripts it pauses; call `report.proceed()` to continue). Nothing is dropped or transformed automatically.

### Step-by-step analysis

```python
import agristatspy as ag
from agristatspy.datasets import load_rcbd_example

df = load_rcbd_example()

print(ag.treatment_summary(df, response="yield", treatment="treatment"))
print(ag.anova_rcbd(df, response="yield", treatment="treatment", block="block"))
print(ag.tukey_hsd(df, response="yield", treatment="treatment"))
```

For a CRD, use `ag.anova_crd(df, response=..., treatment=...)`.

## API overview

| Function | Purpose |
|---|---|
| `ag.design.crd` / `rcbd` / `latin_square` / `factorial` / `split_plot` | Generate a randomized design |
| `ag.treatment_summary(...)` | Per-treatment descriptive statistics and CV% |
| `ag.analyze(...)` | Whole workflow in one call, with export |
| `ag.anova(...)`, `anova_crd`, `anova_rcbd`, `anova_latin_square`, `anova_factorial`, `anova_split_plot` | Analysis of variance |
| `ag.tukey_hsd`, `ag.lsd_test`, `ag.dmrt`, `ag.bonferroni` | Pairwise comparisons (`.letters()` for groupings) |
| `ag.residual_diagnostics(...)` | Residuals, Shapiro–Wilk, Brown–Forsythe |
| `ag.visualization.plot_*` | Layout, means, box, interaction and residual plots |
| `ag.io.write_report(...)`, `ag.io.write_excel(...)` | Markdown / HTML report; multi-sheet Excel workbook |

Design functions share `seed`, `plot_prefix` and `start_plot` (plus design-specific
arguments such as `treatments`, `replications`, `factors`) and return a
`DesignResult` with:

| Member | Description |
|---|---|
| `.table` | Assignment table (a copy) |
| `.metadata` | Immutable design parameters |
| `.validate()` | Check design invariants |
| `.summary()` | Treatment-frequency summary |
| `.text_layout(...)` | Printable field layout |
| `.plot_layout()` | Matplotlib `(fig, ax)` field layout |
| `.to_csv(path)` / `.to_excel(path)` | Export the assignment table |

## Statistical notes

- A significant F-test shows that treatment means are not all equal; it does not say which differ. Follow with a post-hoc comparison.
- Statistical significance does not establish agronomic importance.
- The additive RCBD model does not estimate a treatment-by-block interaction.
- `CV% = √MS_error / grand mean × 100`.
- Observations with missing responses are reported. Blocked and multi-factor analyses require complete, balanced layouts.
- Fisher's LSD does not control the family-wise error rate; diagnostic tests that are not significant do not prove assumptions hold.

Full model descriptions, assumptions, and limitations:
[docs/statistical-methods.md](docs/statistical-methods.md).

## Documentation

- [Installation](docs/installation.md)
- [Quickstart](docs/quickstart.md)
- [Design guide](docs/design-guide.md)
- [Statistical methods](docs/statistical-methods.md)
- [Visualization](docs/visualization.md)
- [Feature comparison](docs/feature-comparison.md)
- [Examples](examples/)

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check src tests examples
```

## Contributing

Contributions are welcome. Please open an issue before submitting a pull
request, especially for new statistical methods. New features should include:

- a defined method, with its assumptions and limitations;
- tests against hand-checkable reference values;
- documentation updates in the same pull request.

## License

Released under the [MIT License](LICENSE).

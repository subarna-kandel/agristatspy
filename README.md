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

- **Design generation** — Completely Randomized Design (CRD) and Randomized Complete Block Design (RCBD) with reproducible, seeded randomization.
- **Field layouts** — printable text layouts with block separators; export to CSV or Excel.
- **Descriptive statistics** — per-treatment summaries and coefficient of variation.
- **Analysis of variance** — one-way ANOVA (CRD) and two-way additive ANOVA (RCBD), with CV% and explicit model assumptions.
- **Mean comparison** — Tukey HSD pairwise comparisons with simultaneous confidence intervals.
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

### Analyze data

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
| `ag.design.crd(...)` | Generate a Completely Randomized Design |
| `ag.design.rcbd(...)` | Generate a Randomized Complete Block Design |
| `ag.treatment_summary(...)` | Per-treatment descriptive statistics and CV% |
| `ag.anova_crd(...)` | One-way ANOVA |
| `ag.anova_rcbd(...)` | Two-way additive ANOVA (treatment + block) |
| `ag.tukey_hsd(...)` | Tukey HSD pairwise comparisons |

Design functions share the parameters `treatments`, `n_treatments`,
`replications`, `seed`, `plot_prefix`, and `start_plot`, and return a
`DesignResult` with:

| Member | Description |
|---|---|
| `.table` | Assignment table (a copy) |
| `.metadata` | Immutable design parameters |
| `.validate()` | Check design invariants |
| `.summary()` | Treatment-frequency summary |
| `.text_layout(...)` | Printable field layout |
| `.to_csv(path)` / `.to_excel(path)` | Export the assignment table |

## Statistical notes

- A significant F-test shows that treatment means are not all equal; it does not say which differ. Follow with a post-hoc comparison.
- Statistical significance does not establish agronomic importance.
- The additive RCBD model does not estimate a treatment-by-block interaction.
- `CV% = √MS_error / grand mean × 100`.
- Observations with missing responses are reported. `anova_rcbd` requires a complete layout (one observation per treatment–block combination).

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

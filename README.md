# AgriStatsPy

**Agricultural experimental design and statistical analysis for Python.**

AgriStatsPy helps agricultural students and researchers create, randomize, visualize, analyze, and report common field-trial designs — without writing lengthy statistical code for every experiment.

---

## Features (Phase 1 — MVP)

| Feature | Status |
|---|---|
| Completely Randomized Design (CRD) | ✅ |
| Randomized Complete Block Design (RCBD) | ✅ |
| Reproducible randomization (seed) | ✅ |
| Text layout with block separators | ✅ |
| CSV export | ✅ |
| Excel export (optional `openpyxl`) | ✅ |
| Descriptive statistics & CV | ✅ |
| One-way ANOVA (CRD) | ✅ |
| Two-way additive ANOVA (RCBD) | ✅ |
| Tukey HSD mean comparisons | ✅ |
| Graphical layout | 🔜 Phase 3 |

---

## Installation

```bash
pip install agristatspy
# With Excel export support:
pip install agristatspy[excel]
```

> **Note:** The package name on PyPI is provisional pending availability checks.

---

## Quick-start

### Design generation

```python
import agristatspy as ag

# Completely Randomized Design
design = ag.design.crd(
    treatments=["Control", "N50", "N100"],
    replications=4,
    seed=12,
)
print(design.text_layout())

# Randomized Complete Block Design
design = ag.design.rcbd(
    treatments=["Control", "N50", "N100", "N150"],
    replications=3,
    seed=42,
)
print(design.text_layout())

# Export the field plan
design.to_csv("field_layout.csv")
```

### Analysis

```python
import pandas as pd
import agristatspy as ag

df = pd.read_csv("observations.csv")

# Descriptive statistics
summary = ag.treatment_summary(df, response="yield", treatment="treatment")
print(summary)

# ANOVA
result = ag.anova_crd(df, response="yield", treatment="treatment")
print(result)

# RCBD ANOVA
result = ag.anova_rcbd(df, response="yield", treatment="treatment", block="block")
print(result)

# Tukey HSD
comparison = ag.tukey_hsd(df, response="yield", treatment="treatment")
print(comparison)
```

### Using built-in example data

```python
from agristatspy.datasets import load_crd_example, load_rcbd_example

df_crd = load_crd_example()
df_rcbd = load_rcbd_example()
```

---

## Design API

### CRD

```python
ag.design.crd(
    treatments=None,     # list of label strings, or None
    n_treatments=None,   # int, used when treatments is not supplied
    replications=1,      # positive int — units per treatment
    seed=None,           # int or None
    plot_prefix="P",     # plot-ID prefix
    start_plot=1,        # first plot number
)
```

### RCBD

```python
ag.design.rcbd(
    treatments=None,     # list of label strings, or None
    n_treatments=None,   # int
    replications=1,      # positive int — number of complete blocks
    seed=None,           # int or None
    plot_prefix="P",     # plot-ID prefix
    start_plot=1,        # first plot number
)
```

### DesignResult object

| Attribute / Method | Description |
|---|---|
| `design.table` | Assignment DataFrame (copy) |
| `design.metadata` | Immutable design parameters |
| `design.validate()` | Run invariant checks |
| `design.summary()` | Human-readable treatment-frequency summary |
| `design.text_layout(...)` | Printable field layout string |
| `design.to_csv(path)` | Export to CSV |
| `design.to_excel(path)` | Export to Excel (requires `openpyxl`) |

---

## Statistical methods

See [`docs/statistical-methods.md`](docs/statistical-methods.md) for full descriptions of the models, assumptions, and limitations.

### Important notices

- A significant F-test shows that treatment means are not all equal; it does not identify which treatments differ. Use Tukey HSD for pairwise comparisons.
- Statistical significance does not establish agronomic importance.
- The standard RCBD model does not estimate a treatment-by-block interaction.
- CV% = √(MS_error) / grand_mean × 100 (agricultural convention).
- Missing observations are reported, not silently removed.

---

## Repository structure

```
agristatspy/
├── src/agristatspy/
│   ├── design/       # CRD, RCBD generators, validation
│   ├── visualization/# Text layout rendering
│   ├── analysis/     # ANOVA, Tukey HSD, descriptive stats
│   ├── io/           # Export helpers
│   ├── datasets/     # Built-in example data
│   └── exceptions.py # Package exception hierarchy
├── tests/            # pytest test suite
├── examples/         # Runnable end-to-end scripts
└── docs/             # Documentation
```

---

## Testing

```bash
pytest
```

---

## License

MIT License — see [LICENSE](LICENSE).

---

## Contributing

Contributions are welcome. Please open an issue before submitting a pull request, especially for new statistical features. All new features must have:

- A defined statistical method, assumptions, and limitations.
- Tests against hand-checkable reference values.
- Documentation updated in the same PR.

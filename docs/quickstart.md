# Quickstart Guide

This quickstart walks you through generating experimental designs and running ANOVA analyses in AgriStatsPy.

---

## 1. Experimental Design

### 1.1 Completely Randomized Design (CRD)

CRD is appropriate when experimental units are homogeneous and no environmental gradient requires blocking.

```python
import agristatspy as ag

# Generate CRD with explicit treatment names
design = ag.design.crd(
    treatments=["Control", "N50", "N100"],
    replications=4,
    seed=12,
)

# View summary & table
print(design.summary())
print(design.table)

# Render formatted text layout
print(design.text_layout())
```

### 1.2 Randomized Complete Block Design (RCBD)

RCBD is used when experimental units can be grouped into homogeneous blocks (e.g., across a slope or soil fertility gradient).

```python
# Generate RCBD with 4 treatments across 3 blocks
rcbd_design = ag.design.rcbd(
    treatments=["Control", "N50", "N100", "N150"],
    replications=3,  # In RCBD, replications equals the number of blocks
    seed=42,
)

# Render formatted layout showing visible block dividers
print(rcbd_design.text_layout(block_separator=True, separator="="))

# Export design to CSV for field crew
rcbd_design.to_csv("rcbd_field_plan.csv")
```

---

## 2. Statistical Analysis

### 2.1 Descriptive Statistics & CV%

```python
from agristatspy.datasets import load_crd_example

data = load_crd_example()

# Summary statistics per treatment
summary = ag.treatment_summary(data, response="yield", treatment="treatment")
print(summary)
```

### 2.2 ANOVA for CRD

```python
# Run one-way ANOVA
result = ag.anova_crd(data, response="yield", treatment="treatment")
print(result)
```

### 2.3 ANOVA for RCBD

```python
from agristatspy.datasets import load_rcbd_example

rcbd_data = load_rcbd_example()

# Run two-way additive ANOVA
rcbd_result = ag.anova_rcbd(
    rcbd_data,
    response="yield",
    treatment="treatment",
    block="block",
)
print(rcbd_result)
```

### 2.4 Post-Hoc Mean Comparisons (Tukey HSD)

```python
# Run Tukey's Honest Significant Difference test
comparisons = ag.tukey_hsd(data, response="yield", treatment="treatment", alpha=0.05)
print(comparisons)
```


---

## 5. More designs and tools

```python
# Latin square, factorial, split-plot
ag.design.latin_square(n_treatments=4, seed=1)
ag.design.factorial({"variety": ["V1", "V2"], "nitrogen": ["N0", "N50"]}, replications=3, layout="rcbd")
ag.design.split_plot(["Rainfed", "Drip"], ["V1", "V2", "V3"], replications=3)

# One entry point for ANOVA
ag.anova(df, "yield", design="rcbd", treatment="treatment", block="block")

# Pairwise comparisons with letters, using the RCBD error term
cmp = ag.tukey_hsd(df, "yield", "treatment", block="block")   # or ag.lsd_test / ag.bonferroni
print(cmp.letters())

# Diagnostics, plots and a report
diag = ag.residual_diagnostics(df, "yield", "treatment", block="block")
fig, ax = ag.visualization.plot_means(df, "yield", "treatment", letters=cmp.letters())
ag.io.write_report("report.md", anova=ag.anova_rcbd(df, "yield", "treatment", "block"), comparison=cmp, diagnostics=diag)
```

See the scripts in `examples/` (05–08) for complete workflows.

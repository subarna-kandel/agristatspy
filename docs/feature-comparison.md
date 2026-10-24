# Feature Comparison: AgriStatsPy vs Other Agricultural & Statistical Tools

This document compares AgriStatsPy's design and analysis workflow with existing tools commonly used in agricultural research, such as R's `agricolae`, `doebioresearch`, and general Python packages like `scipy` and `statsmodels`.

---

## 1. High-Level Comparison

| Capability | AgriStatsPy | R `agricolae` | R `doebioresearch` | Python (`scipy` / `statsmodels`) |
|---|---|---|---|---|
| **Ecosystem** | Python (modern packaging) | R | R | Python |
| **Agricultural Design-to-Analysis Workflow** | Built-in, unified | Built-in | Simplified scripts | Manual data prep & modeling required |
| **Design Generation from Counts** | `ag.design.crd(n_treatments=..., replications=...)` | `design.crd(...)` | N/A (analysis focus) | Manual custom scripts |
| **Reproducible Randomization** | Direct `seed` parameter | `set.seed(...)` | N/A | Manual PRNG handling |
| **Field Text Layout with Separators** | Built-in (`text_layout(block_separator=True)`) | Basic print | N/A | None (user writes ASCII formatter) |
| **Structured Output Object** | Immutable metadata, table, validators | List / data.frame | Text output | Model results wrappers |
| **Tidy Data Integration** | Native Pandas DataFrame | Standard R data.frame | Specific column conventions | DataFrame |
| **ANOVA Model Handling** | Preconfigured CRD / RCBD models | Manual / wrapper | Fixed functions | Explicit formulas (`ols("...").fit()`) |
| **Residual CV% Reporting** | Standard agricultural convention | Calculated | Calculated | Must calculate manually |
| **Mean Comparisons** | Tukey HSD, Fisher LSD, Bonferroni, letter groups | Tukey, LSD, Duncan, SNK | Tukey, LSD, Duncan | `pairwise_tukeyhsd` |
| **Defensive Input Validation** | Strict type & sanity checks | Basic | Basic | Moderate |

---

## 2. Why AgriStatsPy?

1. **Python Native for Agriculture:**
   Agricultural students and researchers increasingly use Python for data science, remote sensing, and GIS, but standard experimental design tools historically lived in R. AgriStatsPy brings a clean, Pythonic experimental design toolchain.

2. **Full Lifecycle (Design → Field Layout → Analysis):**
   Instead of generating plots in one tool and fitting ANOVA models in another, AgriStatsPy supports the entire progression from plot randomization to ANOVA tables and post-hoc tests.

3. **Explicit Field Visualizations:**
   RCBD layouts in field experiments require clear block separation so field staff can see borders and treatment arrangements without ambiguity. AgriStatsPy provides terminal-friendly text layouts with configurable dividers out of the box.

4. **Statistical Rigor & Safety:**
   - Invariant validation ensures no treatment is missing or duplicated within a block.
   - Missing observations are counted and explicitly reported.
   - ANOVA output prints clear model notes emphasizing that statistical significance does not imply agronomic importance.

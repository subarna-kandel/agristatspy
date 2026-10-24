"""
Example 03 — CRD Analysis
===========================
Demonstrates loading a dataset, computing descriptive statistics, running a
one-way ANOVA, and performing Tukey HSD comparisons.
"""
import agristatspy as ag
from agristatspy.datasets import load_crd_example

# -------------------------------------------------------------------
# 1. Load example data
# -------------------------------------------------------------------
df = load_crd_example()
print("=== Data ===")
print(df.to_string(index=False))
print()

# -------------------------------------------------------------------
# 2. Descriptive statistics
# -------------------------------------------------------------------
summary = ag.treatment_summary(df, response="yield", treatment="treatment")
print("=== Descriptive statistics ===")
print(summary.to_string(index=False))
print(f"\nNote: {summary.attrs['cv_note']}\n")

# -------------------------------------------------------------------
# 3. One-way ANOVA
# -------------------------------------------------------------------
result = ag.anova_crd(df, response="yield", treatment="treatment")
print("=== ANOVA ===")
print(result)

# -------------------------------------------------------------------
# 4. Tukey HSD comparisons
# -------------------------------------------------------------------
comparison = ag.tukey_hsd(df, response="yield", treatment="treatment")
print("\n=== Tukey HSD ===")
print(comparison)

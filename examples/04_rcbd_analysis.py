"""
Example 04 — RCBD Analysis
============================
Demonstrates loading RCBD data, running a two-way additive ANOVA (treatment +
block), and performing Tukey HSD comparisons.
"""
import agristatspy as ag
from agristatspy.datasets import load_rcbd_example

# -------------------------------------------------------------------
# 1. Load example data
# -------------------------------------------------------------------
df = load_rcbd_example()
print("=== Data ===")
print(df.to_string(index=False))
print()

# -------------------------------------------------------------------
# 2. Descriptive statistics
# -------------------------------------------------------------------
summary = ag.treatment_summary(df, response="yield", treatment="treatment")
print("=== Descriptive statistics ===")
print(summary.to_string(index=False))
print()

# -------------------------------------------------------------------
# 3. RCBD ANOVA
# -------------------------------------------------------------------
result = ag.anova_rcbd(df, response="yield", treatment="treatment", block="block")
print("=== ANOVA ===")
print(result)

# -------------------------------------------------------------------
# 4. Tukey HSD
# -------------------------------------------------------------------
comparison = ag.tukey_hsd(df, response="yield", treatment="treatment")
print("\n=== Tukey HSD ===")
print(comparison)

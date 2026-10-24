"""
Example 06 — Factorial experiment (variety × nitrogen, RCBD)
============================================================
"""
import agristatspy as ag
from agristatspy.datasets import load_factorial_example

design = ag.design.factorial(
    {"variety": ["V1", "V2"], "nitrogen": ["N0", "N50", "N100"]},
    replications=3,
    layout="rcbd",
    seed=7,
)
print(design.text_layout())
print()

df = load_factorial_example()
print(ag.anova_factorial(df, "yield", ["variety", "nitrogen"], block="block"))

# Compare the six variety × nitrogen combinations, using the RCBD error term.
cmp = ag.tukey_hsd(df, "yield", "treatment", block="block")
print(cmp)
print(cmp.letters())

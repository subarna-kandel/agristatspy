"""
Example 05 — Latin Square
=========================
Generate a randomized 4 × 4 Latin square, then analyze example data.
"""
import agristatspy as ag
from agristatspy.datasets import load_latin_square_example

design = ag.design.latin_square(treatments=["Control", "N50", "N100", "N150"], seed=2026)
print(design.text_layout())
print()

df = load_latin_square_example()
result = ag.anova_latin_square(df, "yield", "treatment", row="row", column="column")
print(result)

# Use the Latin-square error term (rows and columns removed) for comparisons.
print(ag.tukey_hsd(df, "yield", "treatment", block=["row", "column"]))

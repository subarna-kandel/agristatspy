"""
Example 09 — One-command analysis
=================================
Normality check -> ANOVA -> post-hoc letters -> mean ± SEM table, in one call,
then export for a paper or thesis.
"""
import agristatspy as ag
from agristatspy.datasets import load_rcbd_example

df = load_rcbd_example()

report = ag.analyze(df, "yield", "rcbd", treatment="treatment", block="block")
print(report)

# If residuals are not normal, an interactive terminal asks whether to proceed
# or review.  In scripts/notebooks the report is paused instead:
if report.status == "paused":
    report = report.proceed()          # run ANOVA anyway

print(report.thesis_table())           # Mean ± SEM + letters, SEm±, CD, CV%
print(report.anova_table())            # ANOVA table with significance codes

report.export("results.xlsx")          # one sheet per table
report.export("results.tex")           # LaTeX tables
report.export("results.html")          # open in Word

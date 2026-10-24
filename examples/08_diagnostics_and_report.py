"""
Example 08 — Diagnostics, plots, and a written report (RCBD)
============================================================
Writes ``rcbd_report.md``, ``rcbd_report.html`` and PNG figures to the
current directory.
"""
import matplotlib

matplotlib.use("Agg")

import agristatspy as ag  # noqa: E402
from agristatspy.datasets import load_rcbd_example  # noqa: E402

df = load_rcbd_example()

summary = ag.treatment_summary(df, "yield", "treatment")
anova = ag.anova_rcbd(df, "yield", "treatment", "block")
comparison = ag.tukey_hsd(df, "yield", "treatment", block="block")
diagnostics = ag.residual_diagnostics(df, "yield", "treatment", block="block")
print(diagnostics)

fig, _ = ag.visualization.plot_means(
    df, "yield", "treatment", error="se", letters=comparison.letters()
)
fig.savefig("means.png", dpi=200)
fig, _ = ag.visualization.plot_residuals(diagnostics)
fig.savefig("residuals.png", dpi=200)

design = ag.design.rcbd(treatments=["Control", "N50", "N100", "N150"], replications=3, seed=42)
fig, _ = design.plot_layout()
fig.savefig("layout.png", dpi=200)

ag.io.write_excel(
    "rcbd_results.xlsx", design=design, summary=summary, anova=anova,
    comparison=comparison, diagnostics=diagnostics,
)

for ext in ("md", "html"):
    path = ag.io.write_report(
        f"rcbd_report.{ext}",
        title="Nitrogen response trial (RCBD)",
        design=design,
        summary=summary,
        anova=anova,
        comparison=comparison,
        diagnostics=diagnostics,
    )
    print("Wrote", path)

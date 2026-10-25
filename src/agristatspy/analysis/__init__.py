"""
agristatspy.analysis — statistical analysis sub-package.
"""
from agristatspy.analysis.anova import ANOVAResult, anova_crd, anova_rcbd
from agristatspy.analysis.anova_designs import (
    anova,
    anova_factorial,
    anova_latin_square,
    anova_split_plot,
)
from agristatspy.analysis.comparison import (
    ComparisonResult,
    bonferroni,
    dmrt,
    lsd_test,
    tukey_hsd,
)
from agristatspy.analysis.diagnostics import DiagnosticsResult, residual_diagnostics
from agristatspy.analysis.summary import treatment_summary
from agristatspy.analysis.workflow import AnalysisReport, analyze

__all__ = [
    "anova_crd",
    "anova_rcbd",
    "anova",
    "anova_latin_square",
    "anova_factorial",
    "anova_split_plot",
    "tukey_hsd",
    "lsd_test",
    "bonferroni",
    "dmrt",
    "treatment_summary",
    "residual_diagnostics",
    "DiagnosticsResult",
    "analyze",
    "AnalysisReport",
    "ANOVAResult",
    "ComparisonResult",
]

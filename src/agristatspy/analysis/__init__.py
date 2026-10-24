"""
agristatspy.analysis — statistical analysis sub-package.
"""
from agristatspy.analysis.anova import ANOVAResult, anova_crd, anova_rcbd
from agristatspy.analysis.comparison import ComparisonResult, tukey_hsd
from agristatspy.analysis.summary import treatment_summary

__all__ = [
    "anova_crd",
    "anova_rcbd",
    "tukey_hsd",
    "treatment_summary",
    "ANOVAResult",
    "ComparisonResult",
]

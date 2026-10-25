"""
AgriStatsPy — Agricultural experimental design and statistical analysis.

Quick-start
-----------
>>> import agristatspy as ag
>>>
>>> # Generate a CRD
>>> design = ag.design.crd(treatments=["Control", "N50", "N100"], replications=4, seed=12)
>>> print(design.text_layout())
>>>
>>> # Generate an RCBD
>>> design = ag.design.rcbd(
...     treatments=["Control", "N50", "N100", "N150"],
...     replications=3,
...     seed=42,
... )
>>> print(design.text_layout())
>>> design.to_csv("field_layout.csv")
>>>
>>> # Analyze data
>>> import pandas as pd
>>> df = pd.read_csv("observations.csv")
>>> result = ag.analysis.anova_crd(df, response="yield", treatment="treatment")
>>> print(result)
>>> comparison = ag.analysis.tukey_hsd(df, response="yield", treatment="treatment")
>>> print(comparison)
"""
from agristatspy import analysis, design, io, visualization
from agristatspy.analysis import (
    analyze,
    analyze_traits,
    anova,
    anova_alpha_lattice,
    anova_crd,
    anova_factorial,
    anova_latin_square,
    anova_rcbd,
    anova_split_plot,
    bonferroni,
    compare_means,
    correlation,
    dmrt,
    genetic_parameters,
    lsd_test,
    multi_environment,
    pca,
    residual_diagnostics,
    treatment_summary,
    tukey_hsd,
    two_sample_test,
)
from agristatspy.exceptions import (
    AgriStatsError,
    AnalysisError,
    DesignError,
    ExportError,
    ValidationError,
)

__version__ = "0.2.0"
__all__ = [
    "design",
    "analysis",
    "io",
    "visualization",
    "residual_diagnostics",
    "analyze",
    "analyze_traits",
    "anova_alpha_lattice",
    "compare_means",
    "correlation",
    "genetic_parameters",
    "multi_environment",
    "pca",
    "two_sample_test",
    "anova",
    "anova_crd",
    "anova_factorial",
    "anova_latin_square",
    "anova_split_plot",
    "lsd_test",
    "bonferroni",
    "dmrt",
    "anova_rcbd",
    "tukey_hsd",
    "treatment_summary",
    "AgriStatsError",
    "DesignError",
    "ValidationError",
    "AnalysisError",
    "ExportError",
]

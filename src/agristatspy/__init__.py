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
from agristatspy import analysis, design
from agristatspy.analysis import anova_crd, anova_rcbd, treatment_summary, tukey_hsd
from agristatspy.exceptions import (
    AgriStatsError,
    AnalysisError,
    DesignError,
    ExportError,
    ValidationError,
)

__version__ = "0.1.0"
__all__ = [
    "design",
    "analysis",
    "anova_crd",
    "anova_rcbd",
    "tukey_hsd",
    "treatment_summary",
    "AgriStatsError",
    "DesignError",
    "ValidationError",
    "AnalysisError",
    "ExportError",
]

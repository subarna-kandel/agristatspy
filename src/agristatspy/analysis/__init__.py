"""
agristatspy.analysis — statistical analysis sub-package.
"""
from agristatspy.analysis.alpha_lattice import AlphaLatticeResult, anova_alpha_lattice
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
    compare_means,
    dmrt,
    lsd_test,
    tukey_hsd,
)
from agristatspy.analysis.diagnostics import DiagnosticsResult, residual_diagnostics
from agristatspy.analysis.genetics import GeneticParameters, genetic_parameters
from agristatspy.analysis.multi_env import MultiEnvResult, multi_environment
from agristatspy.analysis.multivariate import (
    CorrelationResult,
    PCAResult,
    correlation,
    pca,
)
from agristatspy.analysis.summary import treatment_summary
from agristatspy.analysis.traits import TraitsReport, analyze_traits
from agristatspy.analysis.two_sample import TwoSampleResult, two_sample_test
from agristatspy.analysis.workflow import AnalysisReport, analyze

__all__ = [
    "anova_alpha_lattice",
    "AlphaLatticeResult",
    "multi_environment",
    "MultiEnvResult",
    "analyze_traits",
    "TraitsReport",
    "genetic_parameters",
    "GeneticParameters",
    "correlation",
    "CorrelationResult",
    "pca",
    "PCAResult",
    "two_sample_test",
    "TwoSampleResult",
    "compare_means",
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

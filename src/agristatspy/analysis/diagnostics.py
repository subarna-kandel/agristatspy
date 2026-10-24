"""
Residual diagnostics for CRD and RCBD models.

Diagnostics are screening aids.  A non-significant test does **not** prove that
an assumption holds (especially with small samples), and a significant test
does not by itself justify discarding data.  Outliers are never removed.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from scipy import stats

from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


@dataclass
class DiagnosticsResult:
    """Residuals and assumption-check statistics.

    Attributes
    ----------
    residuals:
        DataFrame with the used rows plus ``fitted`` and ``residual`` columns.
    shapiro_stat, shapiro_p:
        Shapiro–Wilk test of residual normality (``NaN`` if n < 3).
    levene_stat, levene_p:
        Brown–Forsythe (median-centered Levene) test of equal variances
        across treatments (``NaN`` if not computable).
    notes:
        Interpretation caveats.
    """

    residuals: pd.DataFrame
    shapiro_stat: float
    shapiro_p: float
    levene_stat: float
    levene_p: float
    notes: str

    def __str__(self) -> str:
        return "\n".join(
            [
                "Residual diagnostics",
                f"  Shapiro–Wilk (normality):      W = {self.shapiro_stat:.4f}, p = {self.shapiro_p:.4f}",
                f"  Brown–Forsythe (equal var.):   F = {self.levene_stat:.4f}, p = {self.levene_p:.4f}",
                "",
                self.notes,
            ]
        )


def residual_diagnostics(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    block: str | None = None,
) -> DiagnosticsResult:
    """Compute residuals and assumption checks for a CRD or RCBD.

    Parameters
    ----------
    data:
        Tidy DataFrame.
    response, treatment:
        Numeric response and treatment columns.
    block:
        Optional block column; when given, residuals come from the additive
        treatment + block model and the layout must be complete.

    Returns
    -------
    DiagnosticsResult

    Raises
    ------
    AnalysisError
        If columns are invalid, there are fewer than 2 treatments, or a block
        layout is incomplete.
    """
    ctx = "residual_diagnostics"
    cols = [response, treatment] + ([block] if block else [])
    _check_columns(data, cols, ctx)
    _check_numeric(data, response, ctx)
    df = data[cols].dropna(subset=[response]).reset_index(drop=True)
    if df[treatment].nunique() < 2:
        raise AnalysisError(f"{ctx} requires at least 2 treatment groups.")

    y = df[response].astype(float)
    t_mean = df.groupby(treatment)[response].transform("mean").astype(float)
    if block is None:
        fitted = t_mean
    else:
        sizes = df.groupby([treatment, block]).size()
        if (sizes != 1).any() or len(sizes) != df[treatment].nunique() * df[block].nunique():
            raise AnalysisError(
                f"{ctx}: a block layout must have exactly one observation per "
                "treatment-block combination."
            )
        b_mean = df.groupby(block)[response].transform("mean").astype(float)
        fitted = t_mean + b_mean - float(y.mean())

    out = df.copy()
    out["fitted"] = fitted
    out["residual"] = y - fitted

    n = len(out)
    if n >= 3 and out["residual"].nunique() > 1:
        sh = stats.shapiro(out["residual"])
        sh_stat, sh_p = float(sh.statistic), float(sh.pvalue)
    else:
        sh_stat = sh_p = float("nan")

    groups = [g.astype(float).values for _, g in df.groupby(treatment)[response] if len(g) >= 2]
    if len(groups) >= 2:
        lv = stats.levene(*groups, center="median")
        lv_stat, lv_p = float(lv.statistic), float(lv.pvalue)
    else:
        lv_stat = lv_p = float("nan")

    notes = (
        "A non-significant test does not prove the assumption is satisfied,\n"
        "  particularly with few replicates.  Inspect residual plots as well.\n"
        "A significant test is a prompt to investigate (transformation, a\n"
        "  different model, data checks), not a reason to delete observations."
    )
    return DiagnosticsResult(out, sh_stat, sh_p, lv_stat, lv_p, notes)

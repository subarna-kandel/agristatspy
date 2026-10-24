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


def diagnostics_from_fitted(
    df: pd.DataFrame,
    response: str,
    group: str,
    fitted: pd.Series,
) -> DiagnosticsResult:
    """Build a :class:`DiagnosticsResult` from fitted values.

    Parameters
    ----------
    df:
        Rows used in the fit (missing responses already removed).
    response:
        Response column.
    group:
        Column defining the groups for the equal-variance test.
    fitted:
        Fitted values aligned with *df*.
    """
    out = df.copy()
    out["fitted"] = fitted.astype(float).values
    out["residual"] = out[response].astype(float) - out["fitted"]

    if len(out) >= 3 and out["residual"].nunique() > 1:
        sh = stats.shapiro(out["residual"])
        sh_stat, sh_p = float(sh.statistic), float(sh.pvalue)
    else:
        sh_stat = sh_p = float("nan")

    groups = [g.astype(float).values for _, g in out.groupby(group)[response] if len(g) >= 2]
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


def residual_diagnostics(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    block: str | list[str] | None = None,
) -> DiagnosticsResult:
    """Compute residuals and assumption checks for a CRD, RCBD or Latin square.

    Parameters
    ----------
    data:
        Tidy DataFrame.
    response, treatment:
        Numeric response and treatment columns.
    block:
        Optional blocking column, or list of columns (``["row", "column"]``
        for a Latin square).  Residuals then come from the additive
        treatment + blocks model, and the layout must be completely and
        equally crossed.

    Returns
    -------
    DiagnosticsResult

    Raises
    ------
    AnalysisError
        If columns are invalid, there are fewer than 2 treatments, or a
        blocked layout is incomplete.
    """
    ctx = "residual_diagnostics"
    blocks = [] if block is None else ([block] if isinstance(block, str) else list(block))
    cols = [response, treatment, *blocks]
    _check_columns(data, cols, ctx)
    _check_numeric(data, response, ctx)
    df = data[cols].dropna(subset=[response]).reset_index(drop=True)
    if df[treatment].nunique() < 2:
        raise AnalysisError(f"{ctx} requires at least 2 treatment groups.")

    y = df[response].astype(float)
    grand = float(y.mean())
    fitted = pd.Series(grand, index=df.index)
    factors = [treatment, *blocks]
    for i, f in enumerate(factors):
        for g in factors[i + 1:]:
            sizes = df.groupby([f, g]).size()
            if len(sizes) != df[f].nunique() * df[g].nunique() or sizes.nunique() != 1:
                raise AnalysisError(
                    f"{ctx}: '{f}' and '{g}' must be completely and equally crossed "
                    "(e.g. one observation per treatment-block combination)."
                )
    for f in factors:
        fitted = fitted + (df[f].map(df.groupby(f)[response].mean()).astype(float) - grand)
    return diagnostics_from_fitted(df, response, treatment, fitted)

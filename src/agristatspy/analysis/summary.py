"""
Descriptive statistics for agricultural experimental data.

All functions operate on tidy (long-format) DataFrames with one row per
experimental unit.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from agristatspy.exceptions import AnalysisError


def treatment_summary(
    data: pd.DataFrame,
    response: str,
    treatment: str,
) -> pd.DataFrame:
    """Compute per-treatment descriptive statistics.

    Parameters
    ----------
    data:
        Tidy DataFrame with one row per experimental unit.
    response:
        Column name of the numeric response variable.
    treatment:
        Column name of the treatment factor.

    Returns
    -------
    pd.DataFrame
        One row per treatment with columns:
        ``treatment``, ``n``, ``mean``, ``std``, ``se``, ``min``, ``median``,
        ``max``, ``cv_pct``.

        **Coefficient of Variation** (``cv_pct``) is defined as::

            CV% = (residual_SD / grand_mean) × 100

        where ``residual_SD`` is the pooled within-group standard deviation
        (square root of the unweighted mean of within-group variances) and
        ``grand_mean`` is the overall mean of the response.  This matches
        the agricultural convention of reporting ANOVA CV.

        When the grand mean is zero or near-zero (absolute value < 1e-10),
        ``cv_pct`` is set to ``NaN`` and a note is included in the returned
        DataFrame attribute ``attrs["cv_note"]``.

    Raises
    ------
    AnalysisError
        If required columns are missing or the response is not numeric.
    """
    _check_columns(data, [response, treatment], "treatment_summary")
    _check_numeric(data, response, "treatment_summary")

    rows = []
    groups = data.groupby(treatment, sort=True)[response]
    for trt, grp in groups:
        vals = grp.dropna().values.astype(float)
        n = len(vals)
        mean = float(np.mean(vals)) if n > 0 else float("nan")
        std = float(np.std(vals, ddof=1)) if n > 1 else float("nan")
        se = std / math.sqrt(n) if n > 1 else float("nan")
        mn = float(np.min(vals)) if n > 0 else float("nan")
        med = float(np.median(vals)) if n > 0 else float("nan")
        mx = float(np.max(vals)) if n > 0 else float("nan")
        rows.append(
            {
                "treatment": trt,
                "n": n,
                "mean": mean,
                "std": std,
                "se": se,
                "min": mn,
                "median": med,
                "max": mx,
            }
        )

    df = pd.DataFrame(rows)

    # Compute CV using grand mean and pooled within-group SD.
    grand_mean = float(data[response].dropna().mean())
    if abs(grand_mean) < 1e-10:
        df["cv_pct"] = float("nan")
        cv_note = (
            "CV% is undefined because the grand mean is zero or near-zero."
        )
    else:
        # Pooled within-group variance (unweighted mean of within-group variances).
        variances = data.groupby(treatment)[response].var(ddof=1).values
        valid_vars = variances[~np.isnan(variances)]
        pooled_sd = math.sqrt(float(np.mean(valid_vars))) if len(valid_vars) > 0 else float("nan")
        cv_pct = (pooled_sd / abs(grand_mean)) * 100 if not math.isnan(pooled_sd) else float("nan")
        df["cv_pct"] = cv_pct
        cv_note = (
            "CV% = (pooled within-group SD / grand mean) × 100, "
            "consistent with agricultural ANOVA reporting."
        )

    df.attrs["cv_note"] = cv_note
    return df


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _check_columns(data: pd.DataFrame, cols: list[str], context: str) -> None:
    missing = [c for c in cols if c not in data.columns]
    if missing:
        raise AnalysisError(
            f"{context}: required columns missing from data: {missing}. "
            f"Available columns: {list(data.columns)}."
        )


def _check_numeric(data: pd.DataFrame, col: str, context: str) -> None:
    if not pd.api.types.is_numeric_dtype(data[col]):
        raise AnalysisError(
            f"{context}: response column '{col}' must be numeric; "
            f"found dtype {data[col].dtype}."
        )

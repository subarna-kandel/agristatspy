"""
ANOVA for CRD and RCBD designs.

Statistical models
------------------
CRD
    One-way ANOVA (fixed effects):
    y_ij = μ + τ_i + ε_ij,  ε_ij ~ N(0, σ²)

RCBD
    Two-way additive model (treatments + blocks, no interaction):
    y_ij = μ + τ_i + β_j + ε_ij,  ε_ij ~ N(0, σ²)

    The standard additive model does not estimate a treatment-by-block
    interaction.  If interaction is suspected, a more complex design is
    needed.

All ANOVA results are cross-validated against ``scipy.stats`` / ``statsmodels``
calculations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


@dataclass
class ANOVAResult:
    """Container for ANOVA output.

    Attributes
    ----------
    table:
        ANOVA table DataFrame with columns:
        ``source``, ``df``, ``ss``, ``ms``, ``f_stat``, ``p_value``.
    design:
        Short design identifier: ``"CRD"`` or ``"RCBD"``.
    n_obs:
        Total number of observations used.
    grand_mean:
        Grand mean of the response variable.
    cv_pct:
        Coefficient of variation (%) based on the error mean square:
        ``CV% = sqrt(MS_error) / grand_mean × 100``.
        ``NaN`` when the grand mean is zero.
    model_notes:
        Plain-text description of the model and its assumptions.
    missing_reported:
        Number of observations dropped due to missing values.
    """

    table: pd.DataFrame
    design: str
    n_obs: int
    grand_mean: float
    cv_pct: float
    model_notes: str
    missing_reported: int = 0

    def __str__(self) -> str:
        lines = [
            f"ANOVA — {self.design}",
            f"N = {self.n_obs}  |  Grand mean = {self.grand_mean:.4f}"
            f"  |  CV% = {self.cv_pct:.2f}" if not np.isnan(self.cv_pct) else
            f"N = {self.n_obs}  |  Grand mean = {self.grand_mean:.4f}  |  CV% = N/A",
            "",
        ]
        if self.missing_reported:
            lines.insert(1, f"Note: {self.missing_reported} observation(s) with missing response dropped.")
        lines.append(
            self.table.to_string(index=False, float_format=lambda x: f"{x:.4f}")
        )
        lines.append("")
        lines.append("Model notes:")
        lines.append(self.model_notes)
        return "\n".join(lines)


def anova_crd(
    data: pd.DataFrame,
    response: str,
    treatment: str,
) -> ANOVAResult:
    """One-way ANOVA for a Completely Randomized Design.

    Parameters
    ----------
    data:
        Tidy DataFrame with one row per experimental unit.
    response:
        Name of the numeric response column.
    treatment:
        Name of the treatment factor column.

    Returns
    -------
    ANOVAResult

    Raises
    ------
    AnalysisError
        If columns are missing, response is non-numeric, or there are fewer
        than 2 treatment groups or fewer than 2 observations per group.

    Notes
    -----
    The one-way ANOVA F-test assumes:

    - Independent observations.
    - Approximately normally distributed residuals (robust for moderate n).
    - Homogeneity of variance across treatment groups (Levene's test may help
      assess this, though a non-significant test does not prove equality).

    Statistical significance does not establish agronomic importance.
    """
    _check_columns(data, [response, treatment], "anova_crd")
    _check_numeric(data, response, "anova_crd")

    df = data[[treatment, response]].copy()
    n_before = len(df)
    df = df.dropna(subset=[response])
    missing = n_before - len(df)

    groups = df.groupby(treatment)[response]
    group_list = [grp.values.astype(float) for _, grp in groups]
    labels = [name for name, _ in df.groupby(treatment)[response]]

    if len(group_list) < 2:
        raise AnalysisError(
            "anova_crd requires at least 2 treatment groups; "
            f"found {len(group_list)}."
        )
    for lbl, grp in zip(labels, group_list, strict=True):
        if len(grp) < 2:
            raise AnalysisError(
                f"Treatment '{lbl}' has fewer than 2 observations ({len(grp)}); "
                "one-way ANOVA cannot be computed."
            )

    f_stat, p_value = stats.f_oneway(*group_list)

    # Decompose SS manually for the ANOVA table.
    all_vals = df[response].values.astype(float)
    grand_mean = float(np.mean(all_vals))
    n_total = len(all_vals)
    k = len(group_list)

    ss_total = float(np.sum((all_vals - grand_mean) ** 2))
    ss_treatment = float(
        sum(len(g) * (float(np.mean(g)) - grand_mean) ** 2 for g in group_list)
    )
    ss_error = ss_total - ss_treatment

    df_treatment = k - 1
    df_error = n_total - k
    df_total = n_total - 1

    ms_treatment = ss_treatment / df_treatment
    ms_error = ss_error / df_error if df_error > 0 else float("nan")

    cv_pct = (
        (np.sqrt(ms_error) / abs(grand_mean)) * 100
        if abs(grand_mean) > 1e-10
        else float("nan")
    )

    anova_table = pd.DataFrame(
        [
            {
                "source": "Treatment",
                "df": df_treatment,
                "ss": ss_treatment,
                "ms": ms_treatment,
                "f_stat": f_stat,
                "p_value": p_value,
            },
            {
                "source": "Error",
                "df": df_error,
                "ss": ss_error,
                "ms": ms_error,
                "f_stat": float("nan"),
                "p_value": float("nan"),
            },
            {
                "source": "Total",
                "df": df_total,
                "ss": ss_total,
                "ms": float("nan"),
                "f_stat": float("nan"),
                "p_value": float("nan"),
            },
        ]
    )

    notes = (
        "Model: one-way fixed-effects ANOVA (CRD).\n"
        "  y_ij = μ + τ_i + ε_ij,  ε_ij ~ N(0, σ²)\n"
        "Assumptions: independent observations, approximately normal residuals,\n"
        "  and homogeneous variances.\n"
        "A significant F-test indicates treatment means are not all equal;\n"
        "  it does not establish which treatments differ (use a post-hoc test).\n"
        "Statistical significance does not imply agronomic importance.\n"
        "CV% = sqrt(MS_error) / grand_mean × 100."
    )

    return ANOVAResult(
        table=anova_table,
        design="CRD",
        n_obs=n_total,
        grand_mean=grand_mean,
        cv_pct=float(cv_pct),
        model_notes=notes,
        missing_reported=missing,
    )


def anova_rcbd(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    block: str,
) -> ANOVAResult:
    """Two-way additive ANOVA for a Randomized Complete Block Design.

    Parameters
    ----------
    data:
        Tidy DataFrame with one row per experimental unit.
    response:
        Name of the numeric response column.
    treatment:
        Name of the treatment factor column.
    block:
        Name of the block factor column.

    Returns
    -------
    ANOVAResult

    Raises
    ------
    AnalysisError
        If columns are missing, the response is non-numeric, or the data do
        not form a complete block structure (one observation per
        treatment-block combination).

    Notes
    -----
    The standard additive RCBD model is:

        y_ij = μ + τ_i + β_j + ε_ij

    This model does not separately estimate a treatment-by-block interaction.
    If interaction is present, the error term is inflated and F-tests are
    conservative.  Check for interaction by inspecting treatment-by-block
    plots before analysis.

    Statistical significance does not establish agronomic importance.
    """
    _check_columns(data, [response, treatment, block], "anova_rcbd")
    _check_numeric(data, response, "anova_rcbd")

    df = data[[treatment, block, response]].copy()
    n_before = len(df)
    df = df.dropna(subset=[response])
    missing = n_before - len(df)

    # Verify completeness: each (treatment, block) pair must appear exactly once.
    combo_counts = df.groupby([treatment, block]).size()
    if (combo_counts != 1).any():
        bad = combo_counts[combo_counts != 1]
        raise AnalysisError(
            "anova_rcbd expects exactly one observation per treatment-block "
            f"combination.  Problematic combinations:\n{bad}"
        )

    n_total = len(df)
    all_vals = df[response].values.astype(float)
    grand_mean = float(np.mean(all_vals))

    t_labels = df[treatment].unique()
    b_labels = df[block].unique()
    k = len(t_labels)  # number of treatments
    b = len(b_labels)  # number of blocks

    if k < 2:
        raise AnalysisError(
            f"anova_rcbd requires at least 2 treatments; found {k}."
        )
    if b < 2:
        raise AnalysisError(
            f"anova_rcbd requires at least 2 blocks; found {b}."
        )

    # Use statsmodels OLS for the additive model.
    # Rename columns to safe identifiers for the formula.
    safe_df = df.rename(
        columns={response: "response", treatment: "trt", block: "blk"}
    )
    # Ensure factors are treated as categorical.
    safe_df["trt"] = safe_df["trt"].astype(str)
    safe_df["blk"] = safe_df["blk"].astype(str)

    model = ols("response ~ C(trt) + C(blk)", data=safe_df).fit()
    anova_tbl = anova_lm(model, typ=1)

    # Extract SS, df, MS, F, p from statsmodels output.
    # Rows: C(trt), C(blk), Residual
    trt_row = anova_tbl.loc["C(trt)"]
    blk_row = anova_tbl.loc["C(blk)"]
    err_row = anova_tbl.loc["Residual"]

    ss_trt = float(trt_row["sum_sq"])
    df_trt = int(trt_row["df"])
    ms_trt = ss_trt / df_trt
    f_trt = float(trt_row["F"])
    p_trt = float(trt_row["PR(>F)"])

    ss_blk = float(blk_row["sum_sq"])
    df_blk = int(blk_row["df"])
    ms_blk = ss_blk / df_blk
    f_blk = float(blk_row["F"])
    p_blk = float(blk_row["PR(>F)"])

    ss_err = float(err_row["sum_sq"])
    df_err = int(err_row["df"])
    ms_err = ss_err / df_err

    ss_total = float(np.sum((all_vals - grand_mean) ** 2))
    df_total = n_total - 1

    cv_pct = (
        (np.sqrt(ms_err) / abs(grand_mean)) * 100
        if abs(grand_mean) > 1e-10
        else float("nan")
    )

    result_table = pd.DataFrame(
        [
            {
                "source": "Treatment",
                "df": df_trt,
                "ss": ss_trt,
                "ms": ms_trt,
                "f_stat": f_trt,
                "p_value": p_trt,
            },
            {
                "source": "Block",
                "df": df_blk,
                "ss": ss_blk,
                "ms": ms_blk,
                "f_stat": f_blk,
                "p_value": p_blk,
            },
            {
                "source": "Error",
                "df": df_err,
                "ss": ss_err,
                "ms": ms_err,
                "f_stat": float("nan"),
                "p_value": float("nan"),
            },
            {
                "source": "Total",
                "df": df_total,
                "ss": ss_total,
                "ms": float("nan"),
                "f_stat": float("nan"),
                "p_value": float("nan"),
            },
        ]
    )

    notes = (
        "Model: two-way additive fixed-effects ANOVA (RCBD).\n"
        "  y_ij = μ + τ_i + β_j + ε_ij,  ε_ij ~ N(0, σ²)\n"
        "This model does not estimate a treatment-by-block interaction.\n"
        "Assumptions: independent residuals, approximately normal residuals,\n"
        "  additive (no interaction) treatment and block effects.\n"
        "A significant treatment F-test indicates treatment means are not all equal.\n"
        "Statistical significance does not imply agronomic importance.\n"
        "CV% = sqrt(MS_error) / grand_mean × 100."
    )

    return ANOVAResult(
        table=result_table,
        design="RCBD",
        n_obs=n_total,
        grand_mean=grand_mean,
        cv_pct=float(cv_pct),
        model_notes=notes,
        missing_reported=missing,
    )

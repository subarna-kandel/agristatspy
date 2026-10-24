"""
Treatment mean comparisons — Tukey HSD.

Only Tukey HSD is implemented in Phase 1.  Additional methods (LSD, Duncan)
will be added in later phases after their method specifications and limitations
have been reviewed and validated.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


@dataclass
class ComparisonResult:
    """Container for pairwise comparison output.

    Attributes
    ----------
    table:
        DataFrame with columns: ``group1``, ``group2``, ``mean_diff``,
        ``lower_ci``, ``upper_ci``, ``p_adj``, ``reject``.
    method:
        Name of the comparison procedure (e.g. ``"Tukey HSD"``).
    alpha:
        Significance level used.
    model_notes:
        Plain-text description of the method and its assumptions.
    """

    table: pd.DataFrame
    method: str
    alpha: float
    model_notes: str

    def __str__(self) -> str:
        lines = [
            f"Mean comparison — {self.method} (α = {self.alpha})",
            "",
            self.table.to_string(index=False, float_format=lambda x: f"{x:.4f}"),
            "",
            "Notes:",
            self.model_notes,
        ]
        return "\n".join(lines)


def tukey_hsd(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float = 0.05,
) -> ComparisonResult:
    """Tukey Honest Significant Difference (HSD) test.

    Performs all pairwise comparisons of treatment means and adjusts p-values
    to control the family-wise error rate.

    Parameters
    ----------
    data:
        Tidy DataFrame with one row per experimental unit.
    response:
        Name of the numeric response column.
    treatment:
        Name of the treatment factor column.
    alpha:
        Family-wise significance level.  Default is ``0.05``.

    Returns
    -------
    ComparisonResult

    Raises
    ------
    AnalysisError
        If columns are missing, the response is non-numeric, or there are
        fewer than 2 treatment groups.

    Notes
    -----
    Tukey HSD assumes:

    - Balanced or approximately balanced designs (equal group sizes are
      ideal; the Tukey–Kramer extension handles unequal sizes).
    - Independent observations.
    - Approximately normally distributed residuals.
    - Homogeneous variances.

    The implementation delegates to ``statsmodels.stats.multicomp.pairwise_tukeyhsd``,
    which uses the Tukey–Kramer approach for unequal group sizes.

    Pairwise comparisons are the authoritative output.  Letter groupings
    (compact letter displays) can have multiple valid arrangements and are not
    produced in Phase 1.

    Statistical significance does not establish agronomic importance.
    """
    _check_columns(data, [response, treatment], "tukey_hsd")
    _check_numeric(data, response, "tukey_hsd")

    df = data[[treatment, response]].dropna(subset=[response])

    groups = df[treatment].unique()
    if len(groups) < 2:
        raise AnalysisError(
            f"tukey_hsd requires at least 2 treatment groups; found {len(groups)}."
        )

    result = pairwise_tukeyhsd(
        endog=df[response].values.astype(float),
        groups=df[treatment].values,
        alpha=alpha,
    )

    # Build tidy output table.
    summary_data = result.summary().data
    # Row 0 is the header; rows 1+ are data rows.
    header = summary_data[0]
    rows = summary_data[1:]

    # statsmodels columns: group1, group2, meandiff, p-adj, lower, upper, reject
    tbl = pd.DataFrame(rows, columns=header)
    tbl = tbl.rename(
        columns={
            "group1": "group1",
            "group2": "group2",
            "meandiff": "mean_diff",
            "p-adj": "p_adj",
            "lower": "lower_ci",
            "upper": "upper_ci",
            "reject": "reject",
        }
    )
    # Ensure numeric types.
    for col in ["mean_diff", "p_adj", "lower_ci", "upper_ci"]:
        tbl[col] = pd.to_numeric(tbl[col], errors="coerce")

    notes = (
        "Method: Tukey Honest Significant Difference (HSD).\n"
        f"Family-wise error rate is controlled at α = {alpha}.\n"
        "Implementation: statsmodels pairwise_tukeyhsd (Tukey–Kramer for unequal n).\n"
        "Assumptions: independent observations, approximately normal residuals,\n"
        "  homogeneous variances.\n"
        "'reject=True' indicates the null hypothesis of equal means is rejected\n"
        "  at the specified α level.\n"
        "Confidence intervals are simultaneous intervals for pairwise differences.\n"
        "Statistical significance does not establish agronomic importance.\n"
        "Pairwise comparisons are the primary output; compact letter displays\n"
        "  are a convenience summary only."
    )

    return ComparisonResult(
        table=tbl[["group1", "group2", "mean_diff", "lower_ci", "upper_ci", "p_adj", "reject"]],
        method="Tukey HSD",
        alpha=alpha,
        model_notes=notes,
    )

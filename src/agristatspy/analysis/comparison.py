"""
Treatment-mean comparisons — Tukey HSD, Fisher's LSD, Duncan's DMRT, and Bonferroni.

All procedures use the pooled error mean square (MSE) of the experimental
design: the one-way error for a CRD, or the additive treatment + block error
for an RCBD (pass ``block=``).  Pairwise comparisons are the authoritative
output; compact letter displays (see :meth:`ComparisonResult.letters`) are a
convenience summary and can have more than one valid arrangement.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError

_COLUMNS = ["group1", "group2", "mean_diff", "lower_ci", "upper_ci", "p_adj", "reject"]


@dataclass
class ComparisonResult:
    """Container for pairwise comparison output.

    Attributes
    ----------
    table:
        DataFrame with columns: ``group1``, ``group2``, ``mean_diff``
        (mean of ``group2`` minus mean of ``group1``), ``lower_ci``,
        ``upper_ci``, ``p_adj``, ``reject``.
    method:
        Name of the comparison procedure (e.g. ``"Tukey HSD"``).
    alpha:
        Significance level used.
    model_notes:
        Plain-text description of the method and its assumptions.
    means:
        Treatment means (index = treatment label).
    mse:
        Error mean square used by the procedure.
    df_error:
        Error degrees of freedom used by the procedure.
    """

    table: pd.DataFrame
    method: str
    alpha: float
    model_notes: str
    means: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    mse: float = float("nan")
    df_error: int = 0

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

    def letters(self) -> pd.DataFrame:
        """Compact letter display of treatment means.

        Treatments sharing a letter are not significantly different under the
        chosen procedure.  Treatments are ordered from the largest mean to the
        smallest.

        Returns
        -------
        pd.DataFrame
            Columns ``treatment``, ``mean``, ``group``.
        """
        return compact_letters(self.means, self.table)


def compact_letters(means: pd.Series, table: pd.DataFrame) -> pd.DataFrame:
    """Build a compact letter display from pairwise decisions.

    Parameters
    ----------
    means:
        Treatment means indexed by treatment label.
    table:
        Pairwise table with ``group1``, ``group2``, and ``reject`` columns.

    Returns
    -------
    pd.DataFrame
        Columns ``treatment``, ``mean``, ``group``, ordered by decreasing mean.
    """
    order = list(means.sort_values(ascending=False, kind="stable").index)
    pos = {t: i for i, t in enumerate(order)}
    sig: set[tuple[int, int]] = set()
    for g1, g2, rej in zip(table["group1"], table["group2"], table["reject"], strict=True):
        if bool(rej):
            a, b = sorted((pos[g1], pos[g2]))
            sig.add((a, b))

    n = len(order)
    # Maximal runs of consecutive (by mean) treatments with no significant pair.
    runs: list[tuple[int, int]] = []
    for i in range(n):
        j = i
        while j + 1 < n and not any((k, j + 1) in sig for k in range(i, j + 1)):
            j += 1
        runs.append((i, j))
    runs = [
        r for r in runs
        if not any(o != r and o[0] <= r[0] and r[1] <= o[1] for o in runs)
    ]
    runs = sorted(set(runs))

    alphabet = "abcdefghijklmnopqrstuvwxyz"
    labels: list[str] = [
        alphabet[i] if i < len(alphabet) else f"{alphabet[i % 26]}{i // 26}"
        for i in range(len(runs))
    ]
    groups = ["" for _ in order]
    for (lo, hi), lab in zip(runs, labels, strict=True):
        for idx in range(lo, hi + 1):
            groups[idx] += lab

    return pd.DataFrame(
        {
            "treatment": order,
            "mean": [float(means[t]) for t in order],
            "group": groups,
        }
    )


# ---------------------------------------------------------------------------
# Error-term helpers
# ---------------------------------------------------------------------------

def _error_term(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    block: str | list[str] | None,
    context: str,
) -> tuple[pd.Series, pd.Series, float, int]:
    """Return (means, counts, mse, df_error).

    With no blocking factor this is the pooled one-way error.  With one or more
    blocking factors (an RCBD block, or Latin-square row and column) it is the
    error of the additive treatment + blocks model, which requires a balanced,
    orthogonal layout.
    """
    blocks = [] if block is None else ([block] if isinstance(block, str) else list(block))
    cols = [response, treatment, *blocks]
    _check_columns(data, cols, context)
    _check_numeric(data, response, context)

    df = data[cols].dropna(subset=[response])
    grp = df.groupby(treatment)[response]
    means = grp.mean().astype(float)
    counts = grp.size()

    if len(means) < 2:
        raise AnalysisError(
            f"{context} requires at least 2 treatment groups; found {len(means)}."
        )

    y = df[response].astype(float)
    if not blocks:
        resid = y - df[treatment].map(means).astype(float)
        sse = float((resid**2).sum())
        dfe = int(len(df) - len(means))
    else:
        factors = [treatment, *blocks]
        for i, f in enumerate(factors):
            for g in factors[i + 1:]:
                sizes = df.groupby([f, g]).size()
                full = df[f].nunique() * df[g].nunique()
                if len(sizes) != full or sizes.nunique() != 1:
                    raise AnalysisError(
                        f"{context}: '{f}' and '{g}' must be completely and equally "
                        "crossed (e.g. one observation per treatment-block "
                        "combination in a complete RCBD)."
                    )
        grand = float(y.mean())
        fitted = grand
        for f in factors:
            fitted = fitted + (df[f].map(df.groupby(f)[response].mean()).astype(float) - grand)
        sse = float(((y - fitted) ** 2).sum())
        dfe = int(len(df) - 1 - sum(df[f].nunique() - 1 for f in factors))

    if dfe < 1:
        raise AnalysisError(
            f"{context}: no error degrees of freedom; more observations are needed."
        )
    return means, counts, sse / dfe, dfe


def _validate_alpha(alpha: float) -> None:
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
        raise AnalysisError(f"alpha must be a number strictly between 0 and 1; got {alpha!r}.")


def _pairwise(
    method: str,
    means: pd.Series,
    counts: pd.Series,
    mse: float,
    dfe: int,
    alpha: float,
) -> pd.DataFrame:
    k = len(means)
    n_pairs = k * (k - 1) // 2
    rows = []
    rank = {t: i for i, t in enumerate(means.sort_values(ascending=False, kind="stable").index)}
    for g1, g2 in itertools.combinations(sorted(means.index), 2):
        diff = float(means[g2] - means[g1])
        inv_n = 1.0 / counts[g1] + 1.0 / counts[g2]
        if method == "dmrt":
            # Duncan: means r steps apart are tested at 1 - (1 - alpha)^(r - 1).
            r = abs(rank[g1] - rank[g2]) + 1
            se = np.sqrt(mse / 2.0 * inv_n)
            q = abs(diff) / se if se > 0 else np.inf
            sf = float(stats.studentized_range.sf(q, r, dfe)) if np.isfinite(q) else 0.0
            p = 1.0 - (1.0 - sf) ** (1.0 / (r - 1))
            half = float(stats.studentized_range.ppf((1 - alpha) ** (r - 1), r, dfe)) * se
        elif method == "tukey":
            se = np.sqrt(mse / 2.0 * inv_n)
            q = abs(diff) / se if se > 0 else np.inf
            p = float(stats.studentized_range.sf(q, k, dfe)) if np.isfinite(q) else 0.0
            half = float(stats.studentized_range.ppf(1 - alpha, k, dfe)) * se
        else:
            se = np.sqrt(mse * inv_n)
            t = abs(diff) / se if se > 0 else np.inf
            p_raw = float(2 * stats.t.sf(t, dfe)) if np.isfinite(t) else 0.0
            if method == "lsd":
                p = p_raw
                half = float(stats.t.ppf(1 - alpha / 2, dfe)) * se
            else:  # bonferroni
                p = min(1.0, p_raw * n_pairs)
                half = float(stats.t.ppf(1 - alpha / (2 * n_pairs), dfe)) * se
        rows.append(
            {
                "group1": g1,
                "group2": g2,
                "mean_diff": diff,
                "lower_ci": diff - half,
                "upper_ci": diff + half,
                "p_adj": float(min(max(p, 0.0), 1.0)),
                "reject": bool(p < alpha),
            }
        )
    return pd.DataFrame(rows, columns=_COLUMNS)


def _run(
    method: str,
    label: str,
    context: str,
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float,
    block: str | list[str] | None,
    notes_method: str,
) -> ComparisonResult:
    _validate_alpha(alpha)
    means, counts, mse, dfe = _error_term(data, response, treatment, block, context)
    table = _pairwise(method, means, counts, mse, dfe, alpha)
    error_desc = (
        "additive treatment + blocking-factor error (RCBD / Latin square)"
        if block
        else "pooled one-way error (CRD)"
    )
    notes = (
        f"{notes_method}\n"
        f"Error term: {error_desc}; MSE = {mse:.6g}, df = {dfe}.\n"
        "Unequal group sizes are handled with the Kramer adjustment.\n"
        "Assumptions: independent observations, approximately normal residuals,\n"
        "  homogeneous variances.\n"
        "'reject=True' means the null hypothesis of equal means is rejected at\n"
        "  the stated α.  p_adj is the adjusted p-value for this procedure.\n"
        "Statistical significance does not establish agronomic importance.\n"
        "Pairwise comparisons are the primary output; compact letter displays\n"
        "  are a convenience summary only."
    )
    return ComparisonResult(
        table=table,
        method=label,
        alpha=alpha,
        model_notes=notes,
        means=means,
        mse=float(mse),
        df_error=dfe,
    )


# ---------------------------------------------------------------------------
# Public procedures
# ---------------------------------------------------------------------------

def tukey_hsd(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float = 0.05,
    block: str | list[str] | None = None,
) -> ComparisonResult:
    """Tukey Honest Significant Difference (HSD) test.

    Compares all pairs of treatment means while controlling the family-wise
    error rate at *alpha*, using the studentized range distribution (the
    Tukey–Kramer form when group sizes differ).

    Parameters
    ----------
    data:
        Tidy DataFrame with one row per experimental unit.
    response:
        Numeric response column.
    treatment:
        Treatment factor column.
    alpha:
        Family-wise significance level (default ``0.05``).
    block:
        Optional blocking column, or a list of columns (for a Latin square,
        ``["row", "column"]``).  When given, the error term of the additive
        treatment + blocks model is used and the layout must be complete and
        balanced.

    Returns
    -------
    ComparisonResult

    Raises
    ------
    AnalysisError
        If columns are missing, the response is non-numeric, ``alpha`` is
        invalid, there are fewer than 2 groups, or no error degrees of
        freedom are available.

    Notes
    -----
    Statistical significance does not establish agronomic importance.
    """
    return _run(
        "tukey", "Tukey HSD", "tukey_hsd", data, response, treatment, alpha, block,
        "Method: Tukey Honest Significant Difference (HSD).\n"
        f"Family-wise error rate is controlled at α = {alpha}; confidence\n"
        "  intervals are simultaneous for all pairwise differences.",
    )


def lsd_test(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float = 0.05,
    block: str | list[str] | None = None,
) -> ComparisonResult:
    """Fisher's Least Significant Difference (LSD) test.

    Unadjusted pairwise t-tests using the pooled error mean square.

    .. warning::
       LSD does **not** control the family-wise error rate; the chance of at
       least one false positive grows with the number of comparisons.  It is
       conventionally used only after a significant overall F-test, and
       ``p_adj`` here equals the unadjusted p-value.

    Parameters and returns are as for :func:`tukey_hsd`.
    """
    res = _run(
        "lsd", "Fisher LSD", "lsd_test", data, response, treatment, alpha, block,
        "Method: Fisher's Least Significant Difference (unadjusted t-tests).\n"
        "WARNING: the family-wise error rate is NOT controlled.  Use only after a\n"
        "  significant overall F-test and with few planned comparisons.",
    )
    return res


def bonferroni(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float = 0.05,
    block: str | list[str] | None = None,
) -> ComparisonResult:
    """Bonferroni-adjusted pairwise t-tests.

    Each p-value is multiplied by the number of comparisons (capped at 1).
    Conservative, especially with many treatments; Tukey HSD is usually more
    powerful for all-pairs comparisons.

    Parameters and returns are as for :func:`tukey_hsd`.
    """
    return _run(
        "bonferroni", "Bonferroni", "bonferroni", data, response, treatment, alpha, block,
        "Method: Bonferroni-adjusted pairwise t-tests.\n"
        f"Family-wise error rate is controlled at α = {alpha} (conservative).",
    )


def dmrt(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    alpha: float = 0.05,
    block: str | list[str] | None = None,
) -> ComparisonResult:
    """Duncan's Multiple Range Test (DMRT).

    Compares treatment means using the studentized range, with the critical
    value depending on how many ranks apart two means are: means ``r`` steps
    apart are tested at a per-comparison level of ``1 - (1 - alpha)**(r - 1)``.
    Widely used in agricultural reporting.

    .. warning::
       DMRT is more liberal than Tukey HSD and does not control the
       family-wise error rate at *alpha*.  ``p_adj`` is the p-value for the
       rank-specific test, so ``reject`` is ``p_adj < alpha``.

    Parameters and returns are as for :func:`tukey_hsd`.
    """
    return _run(
        "dmrt", "Duncan DMRT", "dmrt", data, response, treatment, alpha, block,
        "Method: Duncan's Multiple Range Test (rank-dependent critical ranges).\n"
        "WARNING: more liberal than Tukey HSD; the family-wise error rate is NOT\n"
        "  controlled at α.  Use after a significant overall F-test.",
    )

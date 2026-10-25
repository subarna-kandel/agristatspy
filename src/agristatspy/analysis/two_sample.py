"""
Two-group comparisons: Welch t-test, Student t-test and Wilcoxon tests.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


@dataclass
class TwoSampleResult:
    """Result of :func:`two_sample_test`.

    ``table`` holds the per-group ``n``, ``mean``, ``sd``, ``sem`` and
    ``median``.  For t-tests ``estimate`` is ``mean(group2) - mean(group1)``
    with its two-sided confidence interval ``ci``; for Wilcoxon tests both are
    ``None``.  ``effect`` is Cohen's d (t-tests) or the rank-biserial
    correlation (Wilcoxon).
    """

    test: str
    groups: tuple[str, str]
    table: pd.DataFrame
    statistic: float
    df: float | None
    p_value: float
    estimate: float | None
    ci: tuple[float, float] | None
    effect: float
    effect_name: str
    alpha: float
    alternative: str
    notes: str

    @property
    def significant(self) -> bool:
        return bool(self.p_value < self.alpha)

    def __str__(self) -> str:
        g1, g2 = self.groups
        lines = [f"{self.test} — {g1} vs {g2} ({self.alternative})", "",
                 self.table.to_string(index=False, float_format=lambda x: f"{x:.4f}"), ""]
        stat = f"statistic = {self.statistic:.4f}" + (f", df = {self.df:.2f}" if self.df else "")
        lines.append(f"{stat}, p = {self.p_value:.4g}  →  "
                     + ("significant" if self.significant else "not significant")
                     + f" at α = {self.alpha}")
        if self.estimate is not None and self.ci is not None:
            lines.append(f"Mean difference ({g2} − {g1}) = {self.estimate:.4f}, "
                         f"{int((1 - self.alpha) * 100)}% CI [{self.ci[0]:.4f}, {self.ci[1]:.4f}]")
        lines.append(f"{self.effect_name} = {self.effect:.3f}")
        lines += ["", self.notes]
        return "\n".join(lines)


def two_sample_test(
    data: pd.DataFrame,
    response: str,
    group: str,
    test: str = "welch",
    pair: str | None = None,
    alpha: float = 0.05,
    alternative: str = "two-sided",
) -> TwoSampleResult:
    """Compare two groups with a t-test or a Wilcoxon test.

    Parameters
    ----------
    data:
        Tidy DataFrame.
    response:
        Numeric response column.
    group:
        Column with exactly two levels.
    test:
        ``"welch"`` (default; unequal variances), ``"student"`` (equal
        variances) or ``"wilcoxon"`` (rank-sum / Mann–Whitney U, or the
        signed-rank test when *pair* is given).
    pair:
        Column identifying matched units (e.g. plot, soil sample ID).  Makes
        the test paired: paired t-test, or Wilcoxon signed-rank.
    alpha:
        Significance level (also sets the confidence level).
    alternative:
        ``"two-sided"``, ``"less"`` or ``"greater"`` (group 2 relative to
        group 1, with groups in sorted order).

    Returns
    -------
    TwoSampleResult

    Raises
    ------
    AnalysisError
        On missing columns, anything other than 2 groups, too few
        observations, or unmatched pairs.

    Notes
    -----
    The Wilcoxon rank-sum test compares distributions, not strictly means or
    medians, unless the two shapes are alike.  For 3 or more groups use
    ANOVA (or Kruskal–Wallis).  Statistical significance does not establish
    agronomic importance.
    """
    if test not in ("welch", "student", "wilcoxon"):
        raise AnalysisError("test must be 'welch', 'student' or 'wilcoxon'.")
    if alternative not in ("two-sided", "less", "greater"):
        raise AnalysisError("alternative must be 'two-sided', 'less' or 'greater'.")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
        raise AnalysisError(f"alpha must be a number strictly between 0 and 1; got {alpha!r}.")
    cols = [response, group] + ([pair] if pair else [])
    _check_columns(data, cols, "two_sample_test")
    _check_numeric(data, response, "two_sample_test")
    df = data[cols].dropna(subset=[response])
    levels = sorted(df[group].dropna().unique(), key=str)
    if len(levels) != 2:
        raise AnalysisError(
            f"two_sample_test needs exactly 2 groups in '{group}'; found {len(levels)}."
        )
    g1, g2 = (str(v) for v in levels)
    a = df[df[group] == levels[0]]
    b = df[df[group] == levels[1]]
    if pair:
        m = a[[pair, response]].merge(b[[pair, response]], on=pair, suffixes=("_1", "_2"))
        if len(m) < 2 or len(m) != len(a) or len(m) != len(b) or a[pair].duplicated().any() \
                or b[pair].duplicated().any():
            raise AnalysisError(
                f"two_sample_test: '{pair}' must identify one matched observation in each group."
            )
        x, y = m[f"{response}_1"].to_numpy(float), m[f"{response}_2"].to_numpy(float)
    else:
        x, y = a[response].to_numpy(float), b[response].to_numpy(float)
        if len(x) < 2 or len(y) < 2:
            raise AnalysisError("two_sample_test needs at least 2 observations per group.")

    tbl = pd.DataFrame(
        {
            "group": [g1, g2], "n": [len(x), len(y)], "mean": [x.mean(), y.mean()],
            "sd": [x.std(ddof=1), y.std(ddof=1)],
            "sem": [x.std(ddof=1) / np.sqrt(len(x)), y.std(ddof=1) / np.sqrt(len(y))],
            "median": [np.median(x), np.median(y)],
        }
    )
    est = ci = dfree = None
    if test == "wilcoxon":
        if pair:
            d = y - x
            if np.all(d == 0):
                raise AnalysisError("two_sample_test: all paired differences are zero.")
            res = stats.wilcoxon(y, x, alternative=alternative)
            name = "Wilcoxon signed-rank test"
            nz = d[d != 0]
            ranks = stats.rankdata(np.abs(nz))
            eff = float((ranks[nz > 0].sum() - ranks[nz < 0].sum()) / ranks.sum())
        else:
            res = stats.mannwhitneyu(y, x, alternative=alternative)
            name = "Wilcoxon rank-sum (Mann–Whitney U) test"
            eff = float(2 * res.statistic / (len(x) * len(y)) - 1)
        stat, p = float(res.statistic), float(res.pvalue)
        effect_name = "Rank-biserial correlation"
        notes = "Rank-based; does not assume normality. Ties are handled by scipy."
    else:
        if pair:
            d = y - x
            res = stats.ttest_rel(y, x, alternative=alternative)
            dfree = float(len(d) - 1)
            se = d.std(ddof=1) / np.sqrt(len(d))
            est = float(d.mean())
            name = "Paired t-test"
            eff = float(d.mean() / d.std(ddof=1)) if d.std(ddof=1) > 0 else float("nan")
            effect_name = "Cohen's d (paired)"
        else:
            eq = test == "student"
            res = stats.ttest_ind(y, x, equal_var=eq, alternative=alternative)
            v1, v2 = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
            se = float(np.sqrt(v1 + v2)) if not eq else float(
                np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1))
                        / (len(x) + len(y) - 2) * (1 / len(x) + 1 / len(y)))
            )
            dfree = float(len(x) + len(y) - 2) if eq else float(
                (v1 + v2) ** 2 / (v1**2 / (len(x) - 1) + v2**2 / (len(y) - 1))
            )
            est = float(y.mean() - x.mean())
            name = "Student t-test" if eq else "Welch two-sample t-test"
            sp = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1))
                         / (len(x) + len(y) - 2))
            eff = float(est / sp) if sp > 0 else float("nan")
            effect_name = "Cohen's d"
        stat, p = float(res.statistic), float(res.pvalue)
        crit = stats.t.ppf(1 - alpha / 2, dfree)
        ci = (est - crit * se, est + crit * se)
        notes = "Assumes approximately normal data" + (
            " and equal variances." if test == "student" and not pair
            else ". Welch does not assume equal variances." if not pair
            else " of the paired differences."
        )
    return TwoSampleResult(
        test=name, groups=(g1, g2), table=tbl, statistic=stat, df=dfree, p_value=p,
        estimate=est, ci=ci, effect=eff, effect_name=effect_name, alpha=alpha,
        alternative=alternative,
        notes=notes + "\nStatistical significance does not establish agronomic importance.",
    )

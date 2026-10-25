"""
Multi-trait analysis: run the one-command workflow over several responses.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from agristatspy.analysis._export import export_tables
from agristatspy.analysis.workflow import AnalysisReport, analyze
from agristatspy.exceptions import AnalysisError


@dataclass
class TraitsReport:
    """Result of :func:`analyze_traits`.

    Attributes
    ----------
    reports:
        ``{trait: AnalysisReport}`` — the full analysis of each trait.
    summary:
        One row per trait: grand mean, CV%, Shapiro–Wilk p and verdict, the
        ANOVA p-value of every effect (``p_<effect>``), whether any effect is
        significant, and the top-ranked treatment with its letter group.
    skipped:
        ``{trait: reason}`` for traits that could not be analysed.
    """

    reports: dict[str, object]
    summary: pd.DataFrame
    skipped: dict[str, str]
    genetics: object | None = None

    def tables(self) -> dict[str, pd.DataFrame]:
        out = {"Summary": self.summary}
        for trait, rep in self.reports.items():
            if isinstance(rep, AnalysisReport):
                if rep.means is None:
                    continue
                out[f"{trait} means"[:31]] = rep.thesis_table()
                out[f"{trait} ANOVA"[:31]] = rep.anova_table()
            else:
                tabs = rep.tables()
                out[f"{trait} means"[:31]] = rep.thesis_table()
                out[f"{trait} ANOVA"[:31]] = next(v for k, v in tabs.items() if "ANOVA" in k)
        if self.genetics is not None:
            out["Genetic parameters"] = self.genetics.table
        return out

    def export(self, path: str) -> object:
        """Export the summary and every trait's tables (``.xlsx``, ``.csv``,
        ``.md``, ``.html``, ``.tex``)."""
        return export_tables(self.tables(), path, "Multi-trait analysis", str(self))

    def __str__(self) -> str:
        lines = [f"Multi-trait analysis — {len(self.reports)} trait(s)", "",
                 self.summary.to_string(index=False, float_format=lambda x: f"{x:.4g}")]
        if self.skipped:
            lines += ["", "Skipped:", *[f"  {k}: {v}" for k, v in self.skipped.items()]]
        if self.genetics is not None:
            lines += ["", str(self.genetics)]
        lines += ["", "Statistical significance does not imply agronomic importance."]
        return "\n".join(lines)


def _top(means: pd.DataFrame, value: str, group: str | None) -> dict:
    top = means.sort_values(value, ascending=False).iloc[0]
    return {"top_treatment": top["treatment"], "top_group": top[group] if group in means else ""}


def _summary_row(trait: str, rep: object) -> dict:
    """One summary line for any result type returned by :func:`analyze`."""
    if isinstance(rep, AnalysisReport):
        row: dict = {
            "trait": trait,
            "grand_mean": rep.stats.get("grand_mean", np.nan),
            "cv_pct": rep.stats.get("cv_pct", np.nan),
            "shapiro_p": rep.normality.shapiro_p,
            "normal": rep.normal,
        }
        pe = rep.stats.get("p_effects", {})
        alpha = rep.stats.get("alpha", 0.05)
        means = rep.means
        top = _top(means, "mean", "group") if means is not None else {}
    else:
        t = rep.anova.table if hasattr(rep, "anova") else rep.combined.table
        gm = rep.anova.grand_mean if hasattr(rep, "anova") else rep.combined.grand_mean
        cv = rep.anova.cv_pct if hasattr(rep, "anova") else rep.combined.cv_pct
        row = {"trait": trait, "grand_mean": gm, "cv_pct": cv}
        if hasattr(rep, "homogeneity"):
            row["bartlett_p"] = rep.homogeneity["p"]
        t = t.set_index("source")
        pe = {s: float(p) for s, p in t["p_value"].items() if pd.notna(p)}
        alpha = rep.alpha
        value = "adj_mean" if hasattr(rep, "anova") else "mean"
        top = _top(rep.means, value, None)
        if rep.comparisons:
            first = next(iter(rep.comparisons.values())).letters().set_index("treatment")["group"]
            top["top_group"] = first.get(top["top_treatment"], "")
    for src, p in pe.items():
        row[f"p_{src}"] = p
    row["significant"] = any(p < alpha for p in pe.values()) if pe else None
    return {**row, **top}


def analyze_traits(
    data: pd.DataFrame,
    responses: list[str] | None,
    design: str,
    on_nonnormal: str = "proceed",
    genetics: bool = False,
    genetic_k: float = 2.06,
    **kwargs: object,
) -> TraitsReport:
    """Run :func:`agristatspy.analyze` on each of several traits.

    Parameters
    ----------
    data:
        Tidy DataFrame.
    responses:
        Trait columns.  ``None`` uses every numeric column that is not used
        as a design factor.
    design:
        Design name, as for :func:`analyze`.
    on_nonnormal:
        As for :func:`analyze`, but defaulting to ``"proceed"`` because a
        batch run cannot pause per trait.  The Shapiro–Wilk verdict of each
        trait is shown in the summary so you can inspect the ones that failed.
    genetics:
        Also compute :func:`genetic_parameters` for all analysed traits
        (design ``"crd"``/``"rcbd"`` only); stored in ``.genetics``.
    genetic_k:
        Selection differential for genetic advance.
    **kwargs:
        Everything else is passed to :func:`analyze` (``treatment``,
        ``block``, ``posthoc``, ``alpha``, …).

    Returns
    -------
    TraitsReport

    Raises
    ------
    AnalysisError
        If no trait can be analysed.  Traits that fail individually (e.g.
        missing values in a blocked design) are listed in ``skipped``.
    """
    if genetics:
        from agristatspy.analysis.workflow import _DESIGNS, _check_genetics_design

        _check_genetics_design(_DESIGNS.get(str(design).lower().replace("-", "_").replace(" ", "_")))
    factor_cols: list[str] = []
    for k in ("treatment", "block", "row", "column", "main", "sub", "replicate"):
        v = kwargs.get(k)
        if isinstance(v, str):
            factor_cols.append(v)
    for k in ("factors", "environment"):
        v = kwargs.get(k) or []
        factor_cols += [v] if isinstance(v, str) else list(v)
    if responses is None:
        responses = [
            c for c in data.columns
            if pd.api.types.is_numeric_dtype(data[c]) and c not in factor_cols
        ]
    if isinstance(responses, str) or not len(responses):
        raise AnalysisError("analyze_traits: 'responses' must be a non-empty list of columns.")

    reports: dict[str, object] = {}
    skipped: dict[str, str] = {}
    rows = []
    for trait in responses:
        try:
            rep = analyze(data, trait, design, on_nonnormal=on_nonnormal, **kwargs)
        except Exception as exc:  # noqa: BLE001 — report per-trait failures, keep going
            skipped[trait] = str(exc)
            continue
        reports[trait] = rep
        row = _summary_row(trait, rep)
        rows.append(row)
    if not rows:
        raise AnalysisError(
            "analyze_traits: no trait could be analysed. "
            + "; ".join(f"{k}: {v}" for k, v in skipped.items())
        )
    gen = None
    if genetics:
        from agristatspy.analysis.workflow import _DESIGNS, _genetics

        key = _DESIGNS.get(str(design).lower().replace("-", "_").replace(" ", "_"))
        gen = _genetics(
            data, list(reports), key, kwargs.get("treatment"), kwargs.get("block"),
            kwargs.get("environment"), genetic_k,
        )
    return TraitsReport(reports=reports, summary=pd.DataFrame(rows), skipped=skipped, genetics=gen)

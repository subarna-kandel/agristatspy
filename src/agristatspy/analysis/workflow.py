"""
One-command analysis workflow.

``analyze()`` runs, in order:

1. a Shapiro–Wilk normality check of the model residuals;
2. if the residuals look normal (or you choose to proceed anyway), the ANOVA
   for the design;
3. if the treatment effect is significant, a post-hoc comparison with
   compact letter groups;
4. summary statistics (means, SD, SEM, SEm±, CD, CV%, grand mean).

If normality is rejected the workflow does not decide for you: it asks
(interactive sessions), or pauses and returns a report you can ``proceed()``
or review.  Nothing is dropped or transformed automatically.
"""
from __future__ import annotations

import html
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis.anova import ANOVAResult
from agristatspy.analysis.anova_designs import anova
from agristatspy.analysis.comparison import (
    ComparisonResult,
    bonferroni,
    dmrt,
    lsd_test,
    tukey_hsd,
)
from agristatspy.analysis.diagnostics import (
    DiagnosticsResult,
    diagnostics_from_fitted,
    residual_diagnostics,
)
from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError, ExportError

_DESIGNS = {
    "crd": "crd", "rcbd": "rcbd", "rbd": "rcbd",
    "latin_square": "latin_square", "latin": "latin_square", "lsq": "latin_square",
    "factorial": "factorial",
    "split_plot": "split_plot", "splitplot": "split_plot",
}
_POSTHOC = {"tukey": tukey_hsd, "lsd": lsd_test, "dmrt": dmrt, "bonferroni": bonferroni}
_POSTHOC_MENU = [*_POSTHOC, "all", "none"]


def _resolve_posthoc(posthoc) -> list[str] | None:
    """Normalise ``posthoc`` to a list of method keys, ``[]`` for none, or
    ``None`` meaning "ask the user".  Raises for unknown names."""
    if posthoc is None:
        return []
    items = [posthoc] if isinstance(posthoc, str) else list(posthoc)
    keys = [str(i).lower().strip() for i in items]
    if keys == ["ask"]:
        return None
    if keys == ["none"]:
        return []
    if "all" in keys:
        return list(_POSTHOC)
    bad = [k for k in keys if k not in _POSTHOC]
    if bad or not keys:
        raise AnalysisError(
            f"posthoc must be one or more of {sorted(_POSTHOC)}, 'all', 'ask' or None; got {posthoc!r}."
        )
    return list(dict.fromkeys(keys))


def _ask_posthoc(ask: Callable[[str], str] | None) -> list[str] | None:
    """Prompt for the post-hoc test(s); ``None`` if no way to ask."""
    if ask is None and sys.stdin is not None and sys.stdin.isatty():
        ask = input
    if ask is None:
        return None
    opts = ", ".join(f"{i}={k}" for i, k in enumerate(_POSTHOC_MENU, 1))
    reply = ask(
        f"Mean comparison [{opts}] (comma-separate to pick several): "
    ).strip().lower()
    chosen: list[str] = []
    for tok in (t.strip() for t in reply.split(",")):
        if tok.isdigit() and 1 <= int(tok) <= len(_POSTHOC_MENU):
            tok = _POSTHOC_MENU[int(tok) - 1]
        chosen.append(tok)
    return _resolve_posthoc(chosen)
_NON_EFFECT = {"Block", "Row", "Column", "Error", "Error (a)", "Error (b)", "Total"}
_CELL_SEP = " | "


# ---------------------------------------------------------------------------
# Report object
# ---------------------------------------------------------------------------

@dataclass
class AnalysisReport:
    """Result of :func:`analyze`.

    Attributes
    ----------
    status:
        ``"complete"`` (ANOVA done), ``"paused"`` (non-normal residuals,
        waiting for a decision), ``"review"`` (you chose to review the data),
        or ``"stopped"``.
    design:
        Design identifier (``"crd"``, ``"rcbd"``, …).
    normality:
        :class:`DiagnosticsResult` (Shapiro–Wilk and Brown–Forsythe).
    normal:
        ``True``/``False`` from the Shapiro–Wilk test at ``normality_alpha``,
        or ``None`` when it could not be computed.
    anova, comparison:
        Results, or ``None`` when not run.  ``comparison`` is the first
        post-hoc selected; ``comparisons`` holds every one run, keyed by
        method (``"tukey"``, ``"lsd"``, ``"dmrt"``, ``"bonferroni"``).
    stats:
        Dict with ``grand_mean``, ``cv_pct``, ``sem_pooled`` (SEm±),
        ``cd`` (critical difference at α), ``alpha``, ``mse``, ``df_error``.
    means:
        Numeric per-treatment table: ``treatment``, ``n``, ``mean``, ``sd``,
        ``sem`` and (if a post-hoc ran) ``group``.
    flagged:
        Observations with large standardized residuals (flagged, never
        removed) when the normality check failed.
    messages:
        Plain-language notes about what was done and why.
    """

    status: str
    design: str
    normality: DiagnosticsResult
    normal: bool | None
    anova: ANOVAResult | None = None
    comparison: ComparisonResult | None = None
    comparisons: dict[str, ComparisonResult] = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    means: pd.DataFrame | None = None
    flagged: pd.DataFrame | None = None
    messages: list[str] = field(default_factory=list)
    decimals: int = 2
    _call: dict = field(default_factory=dict, repr=False)

    # -- control flow -----------------------------------------------------

    def proceed(self) -> AnalysisReport:
        """Continue to ANOVA and post-hoc despite non-normal residuals.

        Returns a new, complete report; this one is unchanged.
        """
        if self.status == "complete":
            return self
        return _execute(**{**self._call, "on_nonnormal": "proceed"})

    # -- tables -----------------------------------------------------------

    def thesis_table(self, decimals: int | None = None, footer: bool = True) -> pd.DataFrame:
        """Publication-style treatment table: ``Mean ± SEM`` plus letters.

        Columns: ``Treatment``, ``Mean ± SEM`` and, if a post-hoc ran,
        ``Group``.  With *footer*, rows for SEm±, CD, CV% and the ANOVA
        p-value are appended.
        """
        self._need_means()
        d = self.decimals if decimals is None else decimals
        m = self.means
        out = pd.DataFrame(
            {
                "Treatment": m["treatment"].astype(str),
                "Mean ± SEM": [
                    f"{mu:.{d}f} ± {se:.{d}f}" if pd.notna(se) else f"{mu:.{d}f}"
                    for mu, se in zip(m["mean"], m["sem"], strict=True)
                ],
            }
        )
        gcols = [c for c in m.columns if c.startswith("group_")]
        if len(gcols) > 1:
            for c in gcols:
                out[f"Group ({self.comparisons[c[6:]].method})"] = m[c].values
        elif "group" in m:
            out["Group"] = m["group"].values
        if footer:
            rows = []
            s = self.stats
            if pd.notna(s.get("sem_pooled", np.nan)):
                rows.append(("SEm ±", f"{s['sem_pooled']:.{d + 1}f}"))
            if pd.notna(s.get("cd", np.nan)):
                rows.append((f"CD ({int(s['alpha'] * 100)}%)", f"{s['cd']:.{d}f}"))
            if pd.notna(s.get("cv_pct", np.nan)):
                rows.append(("CV (%)", f"{s['cv_pct']:.{d}f}"))
            for src, pv in s.get("p_effects", {}).items():
                rows.append((f"ANOVA p ({src})", _fmt_p(pv)))
            for label, val in rows:
                row = {"Treatment": label, "Mean ± SEM": val}
                for c in out.columns[2:]:
                    row[c] = ""
                out = pd.concat([out, pd.DataFrame([row])], ignore_index=True)
        return out

    def anova_table(self, decimals: int | None = None) -> pd.DataFrame:
        """ANOVA table formatted for a paper (with significance codes)."""
        if self.anova is None:
            raise AnalysisError(f"No ANOVA was run (status: {self.status}).")
        d = self.decimals if decimals is None else decimals
        t = self.anova.table
        sig = []
        for _src, p in zip(t["source"], t["p_value"], strict=True):
            sig.append(_stars(p) if pd.notna(p) else "")
        out = pd.DataFrame(
            {
                "Source": t["source"],
                "df": t["df"].astype(int),
                "SS": t["ss"].map(lambda v: f"{v:.{d + 1}f}"),
                "MS": t["ms"].map(lambda v: "" if pd.isna(v) else f"{v:.{d + 1}f}"),
                "F": t["f_stat"].map(lambda v: "" if pd.isna(v) else f"{v:.{d + 1}f}"),
                "p": t["p_value"].map(lambda v: "" if pd.isna(v) else _fmt_p(v)),
                "Sig.": sig,
            }
        )
        return out

    def tables(self) -> dict[str, pd.DataFrame]:
        """All tables as DataFrames, keyed by sheet name."""
        out: dict[str, pd.DataFrame] = {}
        if self.means is not None:
            out["Treatment means"] = self.thesis_table()
            out["Means (numeric)"] = self.means
        if self.anova is not None:
            out["ANOVA"] = self.anova_table()
        if len(self.comparisons) > 1:
            for c in self.comparisons.values():
                out[f"Pairwise {c.method}"[:31]] = c.table
        elif self.comparison is not None:
            out["Pairwise"] = self.comparison.table
        out["Summary"] = pd.DataFrame(
            {
                "statistic": [
                    "Shapiro–Wilk W", "Shapiro–Wilk p", "Brown–Forsythe p",
                    "Grand mean", "CV (%)", "SEm ±", "CD", "alpha",
                ],
                "value": [
                    self.normality.shapiro_stat, self.normality.shapiro_p,
                    self.normality.levene_p, self.stats.get("grand_mean", np.nan),
                    self.stats.get("cv_pct", np.nan), self.stats.get("sem_pooled", np.nan),
                    self.stats.get("cd", np.nan), self.stats.get("alpha", np.nan),
                ],
            }
        )
        return out

    # -- export -----------------------------------------------------------

    def export(self, path: str | Path, fmt: str | None = None) -> Path:
        """Export the tables for a paper or thesis.

        Format from the extension (or *fmt*): ``.xlsx`` (one sheet per table),
        ``.csv`` (the treatment-means table), ``.md``, ``.html`` (opens in
        Word), or ``.tex`` (LaTeX ``tabular`` environments).

        Raises
        ------
        ExportError
            If nothing has been analyzed yet, the format is unsupported, or
            writing fails.
        """
        if self.means is None:
            raise ExportError(
                f"Nothing to export: analysis is '{self.status}'. "
                "Call .proceed() first, or resolve the normality issue."
            )
        out = Path(path).resolve()
        fmt = (fmt or out.suffix.lstrip(".")).lower()
        tabs = self.tables()
        try:
            if fmt == "csv":
                tabs["Treatment means"].to_csv(out, index=False)
            elif fmt in ("xlsx", "excel"):
                try:
                    import openpyxl  # noqa: F401
                except ImportError as exc:
                    raise ExportError(
                        "Excel export requires openpyxl.  Install it with: "
                        "pip install agristatspy[excel]"
                    ) from exc
                from agristatspy.io.excel import _autofit

                with pd.ExcelWriter(out, engine="openpyxl") as w:
                    for name, df in tabs.items():
                        df.to_excel(w, sheet_name=name[:31], index=False)
                    pd.DataFrame({"": str(self).splitlines()}).to_excel(
                        w, sheet_name="Report", index=False, header=False
                    )
                    for ws in w.book.worksheets:
                        _autofit(ws)
            elif fmt in ("md", "markdown"):
                out.write_text(_to_markdown(self, tabs), encoding="utf-8")
            elif fmt == "html":
                out.write_text(_to_html(self, tabs), encoding="utf-8")
            elif fmt in ("tex", "latex"):
                out.write_text(_to_latex(tabs), encoding="utf-8")
            else:
                raise ExportError(
                    f"Unsupported export format {fmt!r}. Use xlsx, csv, md, html or tex."
                )
        except ExportError:
            raise
        except OSError as exc:
            raise ExportError(f"Failed to write {out}: {exc}") from exc
        return out

    # -- display ----------------------------------------------------------

    def _need_means(self) -> None:
        if self.means is None:
            raise AnalysisError(
                f"No results yet (status: {self.status}). Call .proceed() to continue."
            )

    def __str__(self) -> str:
        n = self.normality
        lines = [
            f"AgriStatsPy analysis — {self.design.upper()}   [status: {self.status}]",
            "",
            "1. Normality of residuals",
            f"   Shapiro–Wilk: W = {n.shapiro_stat:.4f}, p = {n.shapiro_p:.4f}  →  "
            + {True: "consistent with normality", False: "NOT normal",
               None: "could not be assessed"}[self.normal],
            f"   Brown–Forsythe (equal variances): p = {n.levene_p:.4f}",
        ]
        if self.messages:
            lines += ["", *[f"   • {m}" for m in self.messages]]
        if self.flagged is not None and len(self.flagged):
            lines += ["", "   Observations with large standardized residuals (kept, not removed):",
                      *["     " + ln for ln in self.flagged.to_string(index=False).splitlines()]]
        if self.anova is not None:
            s = self.stats
            lines += [
                "",
                "2. ANOVA",
                f"   Grand mean = {s['grand_mean']:.4f}   CV = "
                + (f"{s['cv_pct']:.2f}%" if pd.notna(s["cv_pct"]) else "N/A")
                + (f"   SEm± = {s['sem_pooled']:.4f}   CD({int(s['alpha'] * 100)}%) = {s['cd']:.4f}"
                   if pd.notna(s.get("cd", np.nan)) else ""),
                "",
                self.anova_table().to_string(index=False),
            ]
        if self.means is not None:
            names = ", ".join(c.method for c in self.comparisons.values())
            lines += ["", "3. Treatment means (mean ± SEM" + (", letters from "
                      + names if self.comparison else "") + ")", "",
                      self.thesis_table().to_string(index=False)]
            if self.comparison is not None:
                lines += ["", "   Means sharing a letter are not significantly different "
                          f"({names}, α = {self.comparison.alpha})."]
        lines += ["", "Statistical significance does not imply agronomic importance."]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _fmt_p(p: float) -> str:
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def _stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _to_markdown(rep: AnalysisReport, tabs: dict[str, pd.DataFrame]) -> str:
    from agristatspy.io.report import _md_table

    parts = [f"# Analysis — {rep.design.upper()}", ""]
    for name, df in tabs.items():
        parts += [f"## {name}", "", _md_table(df), ""]
    parts += ["## Report", "", "```text", str(rep), "```", ""]
    return "\n".join(parts)


def _to_html(rep: AnalysisReport, tabs: dict[str, pd.DataFrame]) -> str:
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>AgriStatsPy analysis</title>",
        "<style>body{font-family:'Times New Roman',serif;max-width:900px;margin:2em auto}"
        "table{border-collapse:collapse;margin:1em 0}"
        "th,td{padding:4px 10px;text-align:left}"
        "th{border-top:2px solid #000;border-bottom:1px solid #000}"
        "tr:last-child td{border-bottom:2px solid #000}</style></head><body>",
        f"<h1>Analysis — {html.escape(rep.design.upper())}</h1>",
    ]
    for name, df in tabs.items():
        parts.append(f"<h2>{html.escape(name)}</h2>")
        parts.append(df.to_html(index=False, na_rep="", border=0,
                                float_format=lambda x: f"{x:.4f}"))
    parts.append(f"<pre>{html.escape(str(rep))}</pre></body></html>")
    return "\n".join(parts)


def _tex(s: object) -> str:
    out = str(s)
    for ch, rep in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("_", r"\_"),
                    ("#", r"\#"), ("$", r"\$"), ("{", r"\{"), ("}", r"\}"), ("±", r"$\pm$")):
        out = out.replace(ch, rep)
    return out


def _to_latex(tabs: dict[str, pd.DataFrame]) -> str:
    parts = []
    for name, df in tabs.items():
        if name == "Means (numeric)":
            continue
        cols = "l" * len(df.columns)
        body = []
        for row in df.itertuples(index=False):
            cells = [("" if (isinstance(v, float) and pd.isna(v)) else
                      (f"{v:.4f}" if isinstance(v, float) else str(v))) for v in row]
            body.append(" & ".join(_tex(c) for c in cells) + r" \\")
        parts += [
            f"% {name}",
            r"\begin{table}[ht]", r"\centering", f"\\caption{{{_tex(name)}}}",
            f"\\begin{{tabular}}{{{cols}}}", r"\hline",
            " & ".join(_tex(c) for c in df.columns) + r" \\", r"\hline",
            *body, r"\hline", r"\end{tabular}", r"\end{table}", "",
        ]
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Engine
# ---------------------------------------------------------------------------

def _normality(
    key: str, df: pd.DataFrame, response: str, treatment, block, row, column, factors, main, sub
) -> tuple[DiagnosticsResult, pd.DataFrame, str | None]:
    """Return (diagnostics, working data, post-hoc group column)."""
    work = df
    if key == "crd":
        return residual_diagnostics(work, response, treatment), work, treatment
    if key == "rcbd":
        return residual_diagnostics(work, response, treatment, block=block), work, treatment
    if key == "latin_square":
        return residual_diagnostics(work, response, treatment, block=[row, column]), work, treatment
    if key == "factorial":
        work = work.copy()
        work["__cell__"] = work[factors].astype(str).agg(_CELL_SEP.join, axis=1)
        return (
            residual_diagnostics(work, response, "__cell__", block=block),
            work,
            "__cell__",
        )
    # split-plot: fitted = mean(block, A) + mean(A, B) - mean(A)
    work = work.copy()
    for c in (main, sub, block):
        if c not in work.columns:
            raise AnalysisError(f"analyze: column {c!r} not found in data.")
    _check_numeric(work, response, "analyze")
    d = work.dropna(subset=[response]).reset_index(drop=True)
    d["__cell__"] = d[main].astype(str) + _CELL_SEP + d[sub].astype(str)
    cells = d.groupby([block, main, sub]).size()
    n_cells = d[block].nunique() * d[main].nunique() * d[sub].nunique()
    if len(cells) != n_cells or cells.nunique() != 1:
        raise AnalysisError(
            "analyze: a split-plot requires a complete, balanced layout "
            "(one observation per block × main × sub combination)."
        )
    fitted = (
        d.groupby([block, main])[response].transform("mean")
        + d.groupby([main, sub])[response].transform("mean")
        - d.groupby(main)[response].transform("mean")
    )
    return diagnostics_from_fitted(d, response, "__cell__", fitted), d, "__cell__"


def _flag_outliers(diag: DiagnosticsResult, response: str, limit: float = 2.5) -> pd.DataFrame:
    r = diag.residuals
    sd = float(r["residual"].std(ddof=1))
    if not sd or np.isnan(sd):
        return r.iloc[0:0]
    z = r["residual"] / sd
    keep = z.abs() > limit
    out = r.loc[keep].copy()
    out["std_resid"] = z[keep].round(2)
    return out.drop(columns=["fitted"]).reset_index(drop=True)


def _prompt_user(diag: DiagnosticsResult, ask: Callable[[str], str]) -> str:
    q = (
        f"\nResiduals are NOT normal (Shapiro–Wilk p = {diag.shapiro_p:.4f}).\n"
        "  [p] proceed with ANOVA anyway (ANOVA is fairly robust with balanced data)\n"
        "  [r] review the data first (flags large residuals, then stop)\n"
        "Your choice [p/r]: "
    )
    while True:
        ans = str(ask(q)).strip().lower()
        if ans in ("p", "proceed", "y", "yes"):
            return "proceed"
        if ans in ("r", "review", "n", "no"):
            return "review"


def _execute(
    data, response, design, treatment, block, row, column, factors, main, sub,
    alpha, normality_alpha, posthoc, posthoc_when, on_nonnormal, ask, decimals,
) -> AnalysisReport:
    call = {
        "data": data, "response": response, "design": design, "treatment": treatment, "block": block,
        "row": row, "column": column, "factors": factors, "main": main, "sub": sub, "alpha": alpha,
        "normality_alpha": normality_alpha, "posthoc": posthoc, "posthoc_when": posthoc_when,
        "on_nonnormal": on_nonnormal, "ask": ask, "decimals": decimals,
    }
    key = _DESIGNS.get(str(design).lower().replace("-", "_").replace(" ", "_"))
    if key is None:
        raise AnalysisError(
            f"Unknown design {design!r}. Choose from: crd, rcbd, latin_square, factorial, split_plot."
        )
    methods = _resolve_posthoc(posthoc)
    if posthoc_when not in ("significant", "always"):
        raise AnalysisError("posthoc_when must be 'significant' or 'always'.")
    if on_nonnormal not in ("ask", "proceed", "review", "stop"):
        raise AnalysisError("on_nonnormal must be 'ask', 'proceed', 'review' or 'stop'.")
    for name, v in (("alpha", alpha), ("normality_alpha", normality_alpha)):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not 0 < v < 1:
            raise AnalysisError(f"{name} must be a number strictly between 0 and 1; got {v!r}.")

    required = {
        "crd": {"treatment": treatment},
        "rcbd": {"treatment": treatment, "block": block},
        "latin_square": {"treatment": treatment, "row": row, "column": column},
        "factorial": {"factors": factors},
        "split_plot": {"main": main, "sub": sub, "block": block},
    }[key]
    miss = [k for k, v in required.items() if v is None]
    if miss:
        raise AnalysisError(f"analyze(design={design!r}) requires: {', '.join(miss)}.")
    if key == "factorial" and (isinstance(factors, str) or len(factors) < 2):
        raise AnalysisError("analyze: 'factors' must be a list of at least two columns.")
    used = [response] + [
        c for v in (treatment, block, row, column, main, sub) if v for c in [v]
    ] + list(factors or [])
    _check_columns(data, used, "analyze")
    _check_numeric(data, response, "analyze")

    diag, work, group = _normality(key, data, response, treatment, block, row, column, factors, main, sub)
    p = diag.shapiro_p
    normal: bool | None = None if np.isnan(p) else bool(p >= normality_alpha)

    report = AnalysisReport(
        status="paused", design=key, normality=diag, normal=normal,
        decimals=decimals, _call=call,
    )

    if normal is False:
        choice = on_nonnormal
        if choice == "ask":
            if ask is None and sys.stdin is not None and sys.stdin.isatty():
                ask = input
            if ask is not None:
                choice = _prompt_user(diag, ask)
            else:
                report.messages.append(
                    "Residuals are not normal (Shapiro–Wilk p < "
                    f"{normality_alpha}). Analysis paused: call .proceed() to run the "
                    "ANOVA anyway, or review the data (see 'flagged'). "
                    "Use on_nonnormal='proceed' to skip this check."
                )
                report.flagged = _flag_outliers(diag, response)
                return report
        if choice in ("review", "stop"):
            report.status = "review" if choice == "review" else "stopped"
            report.flagged = _flag_outliers(diag, response)
            report.messages.append(
                "Residuals are not normal. Review the data for entry errors or "
                "outliers (flagged below, none removed), consider a transformation "
                "(e.g. log or square root for counts/skewed data) or a "
                "non-parametric test, then re-run. Call .proceed() to continue anyway."
            )
            return report
        report.messages.append(
            "Residuals are not normal; proceeding with ANOVA as requested. "
            "Interpret p-values with caution."
        )
    elif normal is None:
        report.messages.append(
            "Normality could not be assessed (too few observations or no residual variation); proceeding."
        )

    # ---- ANOVA ---------------------------------------------------------
    res = anova(
        data, response, key, treatment=treatment, block=block, row=row,
        column=column, factors=factors, main=main, sub=sub,
    )
    report.anova = res
    report.status = "complete"
    t = res.table.set_index("source")
    effects = [s for s in t.index if s not in _NON_EFFECT]
    # For split-plot only report effect rows; all effect rows have p-values.
    pvals = {s: float(t.loc[s, "p_value"]) for s in effects if pd.notna(t.loc[s, "p_value"])}
    significant = any(v < alpha for v in pvals.values())

    # ---- Post-hoc ------------------------------------------------------
    comp = None
    if key == "split_plot":
        report.messages.append(
            "Post-hoc letters are not produced for split-plot designs: comparisons "
            "need the correct error term for each factor (not yet supported)."
        )
    elif methods == []:
        pass
    elif posthoc_when == "significant" and not significant:
        report.messages.append(
            f"No significant treatment effect at α = {alpha}; post-hoc comparisons were not run."
        )
    else:
        blocks = {"crd": None, "rcbd": block, "latin_square": [row, column], "factorial": block}[key]
        if methods is None:
            methods = _ask_posthoc(ask)
            if methods is None:
                methods = ["tukey"]
                report.messages.append(
                    "posthoc='ask' needs an interactive terminal or an ask= function; "
                    "using Tukey HSD."
                )
        for m in methods:
            report.comparisons[m] = _POSTHOC[m](work, response, group, alpha=alpha, block=blocks)
        comp = report.comparison = next(iter(report.comparisons.values()), None)
        if not significant and {"lsd", "dmrt"} & set(methods):
            report.messages.append(
                "LSD/DMRT run without a significant F-test; interpret with caution."
            )
        if key == "factorial":
            report.messages.append(
                "Letters compare factor-level combinations (cell means). If the interaction "
                "is not significant, main-effect means are usually reported instead."
            )

    # ---- Means and stats ----------------------------------------------
    dfm = work.dropna(subset=[response])
    if key == "split_plot":
        group = "__cell__"
    order = list(dict.fromkeys(dfm[group].tolist()))
    g = dfm.groupby(group)[response]
    means = pd.DataFrame(
        {
            "treatment": order,
            "n": [int(g.size()[k]) for k in order],
            "mean": [float(g.mean()[k]) for k in order],
            "sd": [float(g.std(ddof=1)[k]) if g.size()[k] > 1 else np.nan for k in order],
        }
    )
    means["sem"] = means["sd"] / np.sqrt(means["n"])
    for m, c in report.comparisons.items():
        lt = c.letters()
        lm = dict(zip(lt["treatment"], lt["group"], strict=True))
        means[f"group_{m}"] = means["treatment"].map(lm)
    if comp is not None:
        means["group"] = means[f"group_{next(iter(report.comparisons))}"]
    report.means = means

    s: dict = {"alpha": alpha, "grand_mean": res.grand_mean, "cv_pct": res.cv_pct, "p_effects": pvals}
    if "Error" in t.index:
        mse, dfe = float(t.loc["Error", "ms"]), int(t.loc["Error", "df"])
        r = float(means["n"].mean())
        s.update(
            mse=mse, df_error=dfe, sem_pooled=float(np.sqrt(mse / r)),
            cd=float(stats.t.ppf(1 - alpha / 2, dfe) * np.sqrt(2 * mse / r)),
        )
    report.stats = s
    return report


def analyze(
    data: pd.DataFrame,
    response: str,
    design: str,
    treatment: str | None = None,
    block: str | None = None,
    row: str | None = None,
    column: str | None = None,
    factors: list[str] | None = None,
    main: str | None = None,
    sub: str | None = None,
    alpha: float = 0.05,
    normality_alpha: float = 0.05,
    posthoc: str | list[str] | None = "tukey",
    posthoc_when: str = "significant",
    on_nonnormal: str = "ask",
    ask: Callable[[str], str] | None = None,
    decimals: int = 2,
) -> AnalysisReport:
    """Run the whole analysis with one command.

    Steps: Shapiro–Wilk on the model residuals → (decision if not normal) →
    ANOVA → post-hoc with letter groups if the effect is significant →
    means, SD, SEM, SEm±, CD, CV% and grand mean.

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per experimental unit.
    response:
        Numeric response column.
    design:
        ``"crd"``, ``"rcbd"``, ``"latin_square"``, ``"factorial"`` or
        ``"split_plot"``.
    treatment, block, row, column, factors, main, sub:
        Column names required by the design (as for :func:`agristatspy.anova`).
    alpha:
        Significance level for ANOVA and post-hoc (default 0.05).
    normality_alpha:
        Shapiro–Wilk threshold below which residuals are called non-normal.
    posthoc:
        Mean-comparison test(s): ``"tukey"`` (default), ``"lsd"``, ``"dmrt"``,
        ``"bonferroni"``, a list such as ``["lsd", "dmrt"]``, ``"all"`` (run
        every test and report each side by side), ``"ask"`` (prompt in an
        interactive terminal, or via ``ask=``) or ``None`` to skip.
    posthoc_when:
        ``"significant"`` (default) runs the post-hoc only when an ANOVA
        effect is significant; ``"always"`` runs it regardless.
    on_nonnormal:
        What to do if residuals are not normal: ``"ask"`` (default) prompts
        in an interactive terminal, otherwise pauses the report;
        ``"proceed"`` runs the ANOVA anyway; ``"review"`` stops and flags
        large residuals; ``"stop"`` just stops.
    ask:
        Optional function used instead of ``input`` to get the user's choice
        (receives the question, returns ``"p"`` or ``"r"``).
    decimals:
        Decimal places in the formatted tables.

    Returns
    -------
    AnalysisReport
        Use ``print(report)``, ``report.thesis_table()``,
        ``report.anova_table()``, ``report.proceed()`` and
        ``report.export("results.xlsx")``.

    Raises
    ------
    AnalysisError
        For invalid arguments or data that do not fit the design.

    Notes
    -----
    Nothing is removed or transformed automatically.  A non-significant
    Shapiro–Wilk test does not prove normality, and with few observations the
    test has little power; inspect residual plots as well.  Statistical
    significance does not establish agronomic importance.
    """
    return _execute(
        data=data, response=response, design=design, treatment=treatment, block=block,
        row=row, column=column, factors=factors, main=main, sub=sub, alpha=alpha,
        normality_alpha=normality_alpha, posthoc=posthoc, posthoc_when=posthoc_when,
        on_nonnormal=on_nonnormal, ask=ask, decimals=decimals,
    )

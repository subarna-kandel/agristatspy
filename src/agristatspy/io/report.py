"""
Markdown and HTML report generation.
"""
from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from agristatspy.exceptions import ExportError


def _fmt(v: object) -> str:
    if isinstance(v, float):
        return "" if pd.isna(v) else f"{v:.4f}"
    return str(v)


def _md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(_fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def _blocks(
    title: str,
    design: object | None,
    summary: pd.DataFrame | None,
    anova: object | None,
    comparison: object | None,
    diagnostics: object | None,
) -> list[tuple[str, str, object]]:
    """Return (kind, heading, payload) items; kind is 'table' or 'text'."""
    items: list[tuple[str, str, object]] = []
    if design is not None:
        items.append(("text", "Design", str(design.metadata)))  # type: ignore[attr-defined]
        items.append(("text", "Layout", design.text_layout()))  # type: ignore[attr-defined]
    if summary is not None:
        items.append(("table", "Descriptive statistics", summary))
        note = summary.attrs.get("cv_note")
        if note:
            items.append(("para", "", note))
    if anova is not None:
        a = anova
        head = (
            f"Design: {a.design}. N = {a.n_obs}; grand mean = {a.grand_mean:.4f}; "  # type: ignore[attr-defined]
            f"CV% = {_fmt(float(a.cv_pct))}."  # type: ignore[attr-defined]
        )
        if a.missing_reported:  # type: ignore[attr-defined]
            head += f" {a.missing_reported} observation(s) with missing response dropped."  # type: ignore[attr-defined]
        items.append(("para", "ANOVA", head))
        items.append(("table", "", a.table))  # type: ignore[attr-defined]
        items.append(("text", "Model notes", a.model_notes))  # type: ignore[attr-defined]
    if comparison is not None:
        c = comparison
        items.append(("para", f"Mean comparison — {c.method} (α = {c.alpha})", ""))  # type: ignore[attr-defined]
        items.append(("table", "", c.table))  # type: ignore[attr-defined]
        if len(c.means):  # type: ignore[attr-defined]
            items.append(("table", "Compact letter display", c.letters()))  # type: ignore[attr-defined]
        items.append(("text", "Method notes", c.model_notes))  # type: ignore[attr-defined]
    if diagnostics is not None:
        items.append(("text", "Residual diagnostics", str(diagnostics)))
    return items


def write_report(
    path: str | Path,
    title: str = "AgriStatsPy report",
    design: object | None = None,
    summary: pd.DataFrame | None = None,
    anova: object | None = None,
    comparison: object | None = None,
    diagnostics: object | None = None,
    fmt: str | None = None,
) -> Path:
    """Write a Markdown or HTML report of the supplied results.

    Parameters
    ----------
    path:
        Output file.  The format is taken from the extension (``.md`` or
        ``.html``) unless *fmt* is given.
    title:
        Report title.
    design, summary, anova, comparison, diagnostics:
        Any combination of a ``DesignResult``, a ``treatment_summary`` table,
        an ``ANOVAResult``, a ``ComparisonResult`` and a ``DiagnosticsResult``.
        Sections appear in that order; omitted ones are skipped.
    fmt:
        ``"md"``, ``"html"`` or ``"xlsx"`` (delegates to
        :func:`agristatspy.io.write_excel`).

    Returns
    -------
    Path
        Resolved path of the written file.

    Raises
    ------
    ExportError
        If nothing was supplied, the format is unsupported, or writing fails.
    """
    out = Path(path).resolve()
    fmt = (fmt or out.suffix.lstrip(".") or "md").lower()
    if fmt in ("xlsx", "excel"):
        from agristatspy.io.excel import write_excel

        return write_excel(out, design, summary, anova, comparison, diagnostics)
    if fmt not in ("md", "markdown", "html"):
        raise ExportError(f"Unsupported report format {fmt!r}. Choose 'md', 'html' or 'xlsx'.")
    items = _blocks(title, design, summary, anova, comparison, diagnostics)
    if not items:
        raise ExportError("write_report needs at least one of: design, summary, anova, comparison, diagnostics.")

    if fmt == "html":
        parts = [
            "<!doctype html><html><head><meta charset='utf-8'>",
            f"<title>{html.escape(title)}</title>",
            "<style>body{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em}"
            "table{border-collapse:collapse;margin:1em 0}td,th{border:1px solid #bbb;padding:4px 8px;text-align:right}"
            "pre{background:#f5f5f5;padding:1em;overflow-x:auto}</style></head><body>",
            f"<h1>{html.escape(title)}</h1>",
        ]
        for kind, heading, payload in items:
            if heading:
                parts.append(f"<h2>{html.escape(heading)}</h2>")
            if kind == "table":
                parts.append(payload.to_html(index=False, float_format=lambda x: f"{x:.4f}", na_rep=""))  # type: ignore[attr-defined]
            elif kind == "text":
                parts.append(f"<pre>{html.escape(str(payload))}</pre>")
            elif payload:
                parts.append(f"<p>{html.escape(str(payload))}</p>")
        parts.append("</body></html>")
        content = "\n".join(parts)
    else:
        parts = [f"# {title}", ""]
        for kind, heading, payload in items:
            if heading:
                parts += [f"## {heading}", ""]
            if kind == "table":
                parts += [_md_table(payload), ""]  # type: ignore[arg-type]
            elif kind == "text":
                parts += ["```text", str(payload), "```", ""]
            elif payload:
                parts += [str(payload), ""]
        content = "\n".join(parts)

    try:
        out.write_text(content, encoding="utf-8")
    except OSError as exc:
        raise ExportError(f"Failed to write report to {out}: {exc}") from exc
    return out

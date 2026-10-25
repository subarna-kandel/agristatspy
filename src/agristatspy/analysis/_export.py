"""
Shared table export used by multi-table results (alpha lattice, multi-environment, …).
"""
from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from agristatspy.exceptions import ExportError


def export_tables(
    tables: dict[str, pd.DataFrame], path: str | Path, title: str, text: str = ""
) -> Path:
    """Write *tables* to ``.xlsx`` (one sheet each), ``.csv`` (first table),
    ``.md``, ``.html`` (opens in Word) or ``.tex``."""
    from agristatspy.analysis.workflow import _tex
    from agristatspy.io.report import _md_table

    out = Path(path).resolve()
    fmt = out.suffix.lstrip(".").lower()
    try:
        if fmt == "csv":
            next(iter(tables.values())).to_csv(out, index=False)
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
                for name, df in tables.items():
                    df.to_excel(w, sheet_name=name[:31], index=False)
                if text:
                    pd.DataFrame({"": text.splitlines()}).to_excel(
                        w, sheet_name="Report", index=False, header=False
                    )
                for ws in w.book.worksheets:
                    _autofit(ws)
        elif fmt in ("md", "markdown"):
            parts = [f"# {title}", ""]
            for name, df in tables.items():
                parts += [f"## {name}", "", _md_table(df), ""]
            if text:
                parts += ["## Report", "", "```text", text, "```", ""]
            out.write_text("\n".join(parts), encoding="utf-8")
        elif fmt in ("html", "htm"):
            parts = [
                "<!doctype html><html><head><meta charset='utf-8'>"
                f"<title>{html.escape(title)}</title>"
                "<style>body{font-family:'Times New Roman',serif;max-width:900px;margin:2em auto}"
                "table{border-collapse:collapse;margin:1em 0}th,td{padding:4px 10px;text-align:left}"
                "th{border-top:2px solid #000;border-bottom:1px solid #000}"
                "tr:last-child td{border-bottom:2px solid #000}</style></head><body>",
                f"<h1>{html.escape(title)}</h1>",
            ]
            for name, df in tables.items():
                parts.append(f"<h2>{html.escape(name)}</h2>")
                parts.append(df.to_html(index=False, na_rep="", border=0,
                                        float_format=lambda x: f"{x:.4f}"))
            if text:
                parts.append(f"<pre>{html.escape(text)}</pre>")
            out.write_text("\n".join(parts) + "</body></html>", encoding="utf-8")
        elif fmt in ("tex", "latex"):
            parts = []
            for name, df in tables.items():
                body = [
                    " & ".join("" if (isinstance(v, float) and pd.isna(v)) else
                               (f"{v:.4f}" if isinstance(v, float) else _tex(v)) for v in row)
                    + r" \\"
                    for row in df.itertuples(index=False)
                ]
                parts.append(
                    f"% {name}\n\\begin{{tabular}}{{{'l' * len(df.columns)}}}\n\\hline\n"
                    + " & ".join(_tex(c) for c in df.columns) + " \\\\\n\\hline\n"
                    + "\n".join(body) + "\n\\hline\n\\end{tabular}\n"
                )
            out.write_text("\n".join(parts), encoding="utf-8")
        else:
            raise ExportError(
                f"Unsupported export format {fmt!r}. Use .xlsx, .csv, .md, .html or .tex."
            )
    except ExportError:
        raise
    except Exception as exc:
        raise ExportError(f"Failed to write {out}: {exc}") from exc
    return out

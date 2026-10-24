"""
Excel workbook export for designs and analysis results.

Requires the optional ``openpyxl`` dependency (``pip install agristatspy[excel]``).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from agristatspy.exceptions import ExportError


def _autofit(ws: object) -> None:
    for col in ws.columns:  # type: ignore[attr-defined]
        width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)  # type: ignore[attr-defined]


def _text_sheet(writer: pd.ExcelWriter, name: str, text: str) -> None:
    pd.DataFrame({"": text.splitlines()}).to_excel(writer, sheet_name=name, index=False, header=False)


def write_excel(
    path: str | Path,
    design: object | None = None,
    summary: pd.DataFrame | None = None,
    anova: object | None = None,
    comparison: object | None = None,
    diagnostics: object | None = None,
) -> Path:
    """Write designs and results to a multi-sheet Excel workbook.

    Sheets (only for the objects supplied): ``Design`` and ``Layout``,
    ``Summary``, ``ANOVA`` and ``ANOVA notes``, ``Comparison``, ``Letters``
    and ``Comparison notes``, ``Residuals`` and ``Diagnostics``.  Numeric
    values are written unrounded.

    Parameters
    ----------
    path:
        Destination ``.xlsx`` file.
    design, summary, anova, comparison, diagnostics:
        A ``DesignResult``, a ``treatment_summary`` table, an ``ANOVAResult``,
        a ``ComparisonResult`` and a ``DiagnosticsResult``.

    Returns
    -------
    Path
        Resolved path of the written file.

    Raises
    ------
    ExportError
        If nothing is supplied, ``openpyxl`` is missing, or writing fails.
    """
    if all(x is None for x in (design, summary, anova, comparison, diagnostics)):
        raise ExportError(
            "write_excel needs at least one of: design, summary, anova, comparison, diagnostics."
        )
    try:
        import openpyxl  # noqa: F401
    except ImportError as exc:
        raise ExportError(
            "Excel export requires openpyxl.  Install it with: pip install agristatspy[excel]"
        ) from exc

    out = Path(path).resolve()
    try:
        with pd.ExcelWriter(out, engine="openpyxl") as writer:
            if design is not None:
                design.table.to_excel(writer, sheet_name="Design", index=False)  # type: ignore[attr-defined]
                _text_sheet(writer, "Layout", design.text_layout())  # type: ignore[attr-defined]
            if summary is not None:
                summary.to_excel(writer, sheet_name="Summary", index=False)
            if anova is not None:
                anova.table.to_excel(writer, sheet_name="ANOVA", index=False)  # type: ignore[attr-defined]
                head = (
                    f"Design: {anova.design}\nN = {anova.n_obs}\n"  # type: ignore[attr-defined]
                    f"Grand mean = {anova.grand_mean}\nCV% = {anova.cv_pct}\n"  # type: ignore[attr-defined]
                    f"Missing responses dropped = {anova.missing_reported}\n\n"  # type: ignore[attr-defined]
                )
                _text_sheet(writer, "ANOVA notes", head + anova.model_notes)  # type: ignore[attr-defined]
            if comparison is not None:
                comparison.table.to_excel(writer, sheet_name="Comparison", index=False)  # type: ignore[attr-defined]
                if len(comparison.means):  # type: ignore[attr-defined]
                    comparison.letters().to_excel(writer, sheet_name="Letters", index=False)  # type: ignore[attr-defined]
                _text_sheet(
                    writer,
                    "Comparison notes",
                    f"{comparison.method} (alpha = {comparison.alpha})\n\n{comparison.model_notes}",  # type: ignore[attr-defined]
                )
            if diagnostics is not None:
                diagnostics.residuals.to_excel(writer, sheet_name="Residuals", index=False)  # type: ignore[attr-defined]
                _text_sheet(writer, "Diagnostics", str(diagnostics))
            for ws in writer.book.worksheets:
                _autofit(ws)
    except ExportError:
        raise
    except Exception as exc:
        raise ExportError(f"Failed to write Excel workbook to {out}: {exc}") from exc
    return out

"""
Export helpers (thin wrappers around DesignResult methods).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from agristatspy.exceptions import ExportError


def export_design(design: object, path: str | Path, fmt: str = "csv", **kwargs: Any) -> Path:
    """Export a design to the requested format.

    Parameters
    ----------
    design:
        A :class:`~agristatspy.design.base.DesignResult` instance.
    path:
        Destination file path.
    fmt:
        ``"csv"`` (default) or ``"excel"``.
    **kwargs:
        Passed to the underlying export method.

    Returns
    -------
    Path
        Resolved path of the written file.
    """
    if fmt == "csv":
        return design.to_csv(path, **kwargs)  # type: ignore[union-attr]
    elif fmt in ("excel", "xlsx"):
        return design.to_excel(path, **kwargs)  # type: ignore[union-attr]
    else:
        raise ExportError(f"Unsupported export format {fmt!r}. Choose 'csv' or 'excel'.")

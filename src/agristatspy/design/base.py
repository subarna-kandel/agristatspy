"""
Base design object returned by all design-generation functions.

The DesignResult class is the stable public interface that callers interact
with after generating a CRD or RCBD.  It intentionally stores the assignment
table as a read-only-by-convention copy so that display operations cannot
silently mutate the design.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from agristatspy.exceptions import ExportError

if TYPE_CHECKING:
    pass


@dataclass(frozen=True)
class DesignMetadata:
    """Immutable summary of design parameters.

    Attributes
    ----------
    design_type:
        Short identifier such as ``"CRD"`` or ``"RCBD"``.
    treatments:
        Ordered list of unique treatment labels.
    n_treatments:
        Number of distinct treatments.
    replications:
        For CRD: replicate units per treatment.
        For RCBD: number of complete blocks.
    total_plots:
        Total number of experimental units.
    seed:
        Random seed used, or ``None`` if no seed was set.
    plot_prefix:
        String prefix used for plot-ID labels.
    start_plot:
        Integer offset for the first plot ID.
    extra:
        Any design-specific extra metadata.
    """

    design_type: str
    treatments: tuple[str, ...]
    n_treatments: int
    replications: int
    total_plots: int
    seed: int | None
    plot_prefix: str
    start_plot: int
    extra: dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [
            f"Design      : {self.design_type}",
            f"Treatments  : {', '.join(self.treatments)}",
            f"Replications: {self.replications}",
            f"Total plots : {self.total_plots}",
            f"Seed        : {self.seed}",
        ]
        for k, v in self.extra.items():
            lines.append(f"{k:<12}: {v}")
        return "\n".join(lines)


class DesignResult:
    """Container returned by :func:`crd` and :func:`rcbd`.

    Parameters
    ----------
    table:
        Assignment DataFrame.  A copy is stored so external mutations do not
        affect the design object.
    metadata:
        Immutable :class:`DesignMetadata` record.
    """

    def __init__(self, table: pd.DataFrame, metadata: DesignMetadata) -> None:
        # Store a copy to prevent external mutation from affecting the design.
        self._table: pd.DataFrame = table.copy()
        self._metadata: DesignMetadata = metadata

    # ------------------------------------------------------------------
    # Read-only properties
    # ------------------------------------------------------------------

    @property
    def table(self) -> pd.DataFrame:
        """Return a copy of the assignment DataFrame."""
        return self._table.copy()

    @property
    def metadata(self) -> DesignMetadata:
        """Return the immutable design metadata."""
        return self._metadata

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> bool:
        """Run design-constraint checks.

        Returns
        -------
        bool
            ``True`` when all checks pass.

        Raises
        ------
        ValidationError
            With a descriptive message when any check fails.
        """
        from agristatspy.design._validators import validate_design

        validate_design(self)
        return True

    # ------------------------------------------------------------------
    # Summaries
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a brief human-readable summary of treatments and counts."""
        lines = [str(self._metadata), ""]
        counts = (
            self._table.groupby("treatment", sort=False)
            .size()
            .rename("count")
            .reset_index()
        )
        lines.append("Treatment frequencies:")
        for _, row in counts.iterrows():
            lines.append(f"  {row['treatment']:<20} {row['count']}")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Text layout
    # ------------------------------------------------------------------

    def text_layout(
        self,
        mode: str = "auto",
        row_length: int | None = None,
        block_separator: bool = True,
        separator: str = "=",
        cell_width: int = 12,
    ) -> str:
        """Return a printable text layout of the field design.

        Parameters
        ----------
        mode:
            ``"auto"`` chooses ``"blocks"`` for RCBD and ``"grid"`` for CRD.
            Other valid values: ``"blocks"``, ``"grid"``, ``"sequence"``.
        row_length:
            Plots per row in ``"grid"`` mode.  Defaults to the number of
            treatments for RCBD blocks.
        block_separator:
            When ``True`` (default) a prominent divider is printed between
            blocks in RCBD mode.
        separator:
            Character used for horizontal rules and block-separator lines.
        cell_width:
            Minimum character width of each cell in a grid/block row.

        Returns
        -------
        str
            Printable layout string.
        """
        from agristatspy.visualization.text import render_layout

        return render_layout(
            design=self,
            mode=mode,
            row_length=row_length,
            block_separator=block_separator,
            separator=separator,
            cell_width=cell_width,
        )

    def __repr__(self) -> str:
        m = self._metadata
        return (
            f"<DesignResult type={m.design_type!r} "
            f"treatments={m.n_treatments} replications={m.replications} "
            f"plots={m.total_plots}>"
        )

    def __str__(self) -> str:
        return self.text_layout()

    # ------------------------------------------------------------------
    # Graphical layout (Phase 3 — stub in Phase 1)
    # ------------------------------------------------------------------

    def plot_layout(self, **kwargs: Any) -> tuple:
        """Return a Matplotlib (figure, axes) for the field layout.

        .. note::
           Full graphical output is implemented in Phase 3.  This stub raises
           ``NotImplementedError`` with a helpful message.
        """
        raise NotImplementedError(
            "Graphical layout (plot_layout) will be available in Phase 3. "
            "Use text_layout() for a terminal-friendly representation."
        )

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def to_csv(self, path: str | Path, **kwargs: Any) -> Path:
        """Export the assignment table to a CSV file.

        Parameters
        ----------
        path:
            Destination file path.
        **kwargs:
            Additional keyword arguments passed to :func:`pandas.DataFrame.to_csv`.

        Returns
        -------
        Path
            Resolved path of the written file.

        Raises
        ------
        ExportError
            If writing fails.
        """
        resolved = Path(path).resolve()
        try:
            self._table.to_csv(resolved, index=False, **kwargs)
        except Exception as exc:
            raise ExportError(f"Failed to write CSV to {resolved}: {exc}") from exc
        return resolved

    def to_excel(self, path: str | Path, sheet_name: str = "Design", **kwargs: Any) -> Path:
        """Export the assignment table to an Excel file.

        Requires the optional ``openpyxl`` dependency::

            pip install agristatspy[excel]

        Parameters
        ----------
        path:
            Destination file path.
        sheet_name:
            Worksheet name.
        **kwargs:
            Additional keyword arguments passed to
            :func:`pandas.DataFrame.to_excel`.

        Returns
        -------
        Path
            Resolved path of the written file.

        Raises
        ------
        ExportError
            If ``openpyxl`` is not installed or writing fails.
        """
        try:
            import openpyxl  # noqa: F401  — check availability
        except ImportError as exc:
            raise ExportError(
                "Excel export requires openpyxl.  Install it with: "
                "pip install agristatspy[excel]"
            ) from exc

        resolved = Path(path).resolve()
        try:
            self._table.to_excel(
                resolved, index=False, sheet_name=sheet_name, **kwargs
            )
        except Exception as exc:
            raise ExportError(
                f"Failed to write Excel file to {resolved}: {exc}"
            ) from exc
        return resolved

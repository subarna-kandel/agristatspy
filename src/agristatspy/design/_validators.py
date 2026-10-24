"""
Internal design-invariant validators.

These are called by :meth:`DesignResult.validate` and optionally during
generation.  They raise :class:`ValidationError` on failure.
"""
from __future__ import annotations

from agristatspy.exceptions import ValidationError


def validate_design(design: object) -> None:
    """Dispatch to the correct validator based on ``design.metadata.design_type``."""
    dtype = design.metadata.design_type  # type: ignore[attr-defined]
    if dtype == "CRD":
        _validate_crd(design)
    elif dtype == "RCBD":
        _validate_rcbd(design)
    else:
        raise ValidationError(f"No validator defined for design type {dtype!r}.")


def _validate_crd(design: object) -> None:
    """Validate CRD invariants."""
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]

    expected_total = meta.n_treatments * meta.replications
    if len(tbl) != expected_total:
        raise ValidationError(
            f"CRD expected {expected_total} experimental units "
            f"({meta.n_treatments} × {meta.replications}), "
            f"but found {len(tbl)}."
        )

    if tbl["plot_id"].nunique() != len(tbl):
        raise ValidationError("Plot IDs are not unique in the CRD design.")

    counts = tbl.groupby("treatment").size()
    for treatment in meta.treatments:
        n = counts.get(treatment, 0)
        if n != meta.replications:
            raise ValidationError(
                f"Treatment '{treatment}' appears {n} times; "
                f"expected {meta.replications}."
            )


def _validate_rcbd(design: object) -> None:
    """Validate RCBD invariants."""
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    n_blocks = meta.replications  # replications == number of blocks in RCBD

    expected_total = meta.n_treatments * n_blocks
    if len(tbl) != expected_total:
        raise ValidationError(
            f"RCBD expected {expected_total} plots "
            f"({meta.n_treatments} treatments × {n_blocks} blocks), "
            f"but found {len(tbl)}."
        )

    if tbl["plot_id"].nunique() != len(tbl):
        raise ValidationError("Plot IDs are not unique in the RCBD design.")

    blocks = tbl["block"].unique()
    if len(blocks) != n_blocks:
        raise ValidationError(
            f"Expected {n_blocks} blocks but found {len(blocks)}."
        )

    for block_id in blocks:
        block_df = tbl[tbl["block"] == block_id]
        block_treatments = set(block_df["treatment"].tolist())
        expected_treatments = set(meta.treatments)
        if block_treatments != expected_treatments:
            missing = expected_treatments - block_treatments
            extra = block_treatments - expected_treatments
            msg = f"Block '{block_id}' does not contain all treatments."
            if missing:
                msg += f" Missing: {sorted(missing)}."
            if extra:
                msg += f" Unexpected: {sorted(extra)}."
            raise ValidationError(msg)

        if len(block_df) != meta.n_treatments:
            raise ValidationError(
                f"Block '{block_id}' has {len(block_df)} plots; "
                f"expected {meta.n_treatments}."
            )

"""
Internal design-invariant validators.

These are called by :meth:`DesignResult.validate` and optionally during
generation.  They raise :class:`ValidationError` on failure.
"""
from __future__ import annotations

from agristatspy.exceptions import ValidationError


def validate_design(design: object) -> None:
    """Dispatch to the correct validator based on ``design.metadata.design_type``."""
    meta = design.metadata  # type: ignore[attr-defined]
    dtype = meta.extra.get("base_design", meta.design_type)
    if dtype == "CRD":
        _validate_crd(design)
    elif dtype == "RCBD":
        _validate_rcbd(design)
    elif dtype == "LSQ":
        _validate_latin_square(design)
    elif dtype == "SplitPlot":
        _validate_split_plot(design)
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


def _validate_latin_square(design: object) -> None:
    """Validate Latin square invariants."""
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    t = meta.n_treatments

    if len(tbl) != t * t:
        raise ValidationError(f"Latin square expected {t * t} plots; found {len(tbl)}.")
    if tbl["plot_id"].nunique() != len(tbl):
        raise ValidationError("Plot IDs are not unique in the Latin square design.")

    expected = set(meta.treatments)
    for axis in ("row", "column"):
        for key, grp in tbl.groupby(axis):
            if len(grp) != t or set(grp["treatment"]) != expected:
                raise ValidationError(
                    f"{axis.capitalize()} {key} does not contain each treatment exactly once."
                )


def _validate_split_plot(design: object) -> None:
    """Validate split-plot invariants."""
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    main_name = meta.extra["main_name"]
    sub_name = meta.extra["sub_name"]
    main = set(meta.extra["main_levels"])
    sub = set(meta.extra["sub_levels"])
    n_blocks = meta.replications

    if len(tbl) != len(main) * len(sub) * n_blocks:
        raise ValidationError(
            f"Split-plot expected {len(main) * len(sub) * n_blocks} plots; found {len(tbl)}."
        )
    if tbl["plot_id"].nunique() != len(tbl):
        raise ValidationError("Plot IDs are not unique in the split-plot design.")

    for block_id, bdf in tbl.groupby("block"):
        mains = bdf.groupby("main_plot")[main_name].nunique()
        if (mains != 1).any():
            raise ValidationError(f"A main plot in block {block_id} has mixed main-plot levels.")
        if set(bdf[main_name]) != main or bdf["main_plot"].nunique() != len(main):
            raise ValidationError(
                f"Block {block_id} does not contain each main-plot level exactly once."
            )
    for mp, mdf in tbl.groupby("main_plot"):
        if len(mdf) != len(sub) or set(mdf[sub_name]) != sub:
            raise ValidationError(
                f"Main plot {mp} does not contain each sub-plot level exactly once."
            )

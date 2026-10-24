"""
Text-layout renderer for CRD and RCBD designs.

All rendering functions accept a :class:`~agristatspy.design.base.DesignResult`
and return a plain string suitable for printing, testing, or saving to a file.
Rendering never alters the underlying assignment table.
"""
from __future__ import annotations


def render_layout(
    design: object,
    mode: str = "auto",
    row_length: int | None = None,
    block_separator: bool = True,
    separator: str = "=",
    cell_width: int = 12,
) -> str:
    """Dispatch to the correct renderer based on design type and mode.

    Parameters
    ----------
    design:
        A :class:`~agristatspy.design.base.DesignResult` instance.
    mode:
        ``"auto"`` selects ``"blocks"`` for RCBD and ``"grid"`` for CRD.
        Explicit values: ``"blocks"``, ``"grid"``, ``"sequence"``.
    row_length:
        Plots per row in grid/blocks mode.
    block_separator:
        Print a prominent separator between blocks (RCBD).
    separator:
        Character for horizontal rules.
    cell_width:
        Minimum cell width.

    Returns
    -------
    str
        Formatted layout string.
    """
    meta = design.metadata  # type: ignore[attr-defined]
    dtype = meta.extra.get("base_design", meta.design_type)
    resolved_mode = mode

    if mode == "table":
        return design.table.to_string(index=False)  # type: ignore[attr-defined]

    if row_length is not None and (
        isinstance(row_length, bool) or not isinstance(row_length, int) or row_length < 1
    ):
        raise ValueError(f"row_length must be a positive integer; got {row_length!r}.")

    if mode == "auto":
        resolved_mode = "blocks" if dtype in ("RCBD", "SplitPlot") else "grid"

    if dtype == "CRD":
        return _render_crd(
            design,
            mode=resolved_mode,
            row_length=row_length,
            separator=separator,
            cell_width=cell_width,
        )
    elif dtype == "RCBD":
        return _render_rcbd(
            design,
            mode=resolved_mode,
            row_length=row_length,
            block_separator=block_separator,
            separator=separator,
            cell_width=cell_width,
        )
    elif dtype == "LSQ":
        return _render_latin_square(design, cell_width=cell_width)
    elif dtype == "SplitPlot":
        return _render_split_plot(
            design,
            block_separator=block_separator,
            separator=separator,
            cell_width=cell_width,
        )
    else:
        raise ValueError(f"No text renderer for design type {dtype!r}.")


# ---------------------------------------------------------------------------
# Shared cell/row helpers
# ---------------------------------------------------------------------------

def _cell(content: str, width: int) -> str:
    """Format *content* into a padded cell (without leading pipe)."""
    return f" {content:<{width}} "


def _hrule(n_cols: int, cell_width: int, sep: str = "-") -> str:
    """Horizontal rule spanning *n_cols* columns."""
    col = sep * (cell_width + 2)
    return "+" + "+".join([col] * n_cols) + "+"


def _row_line(values: list[str], cell_width: int) -> str:
    """Build one data row with pipe separators."""
    cells = [_cell(v, cell_width) for v in values]
    return "|" + "|".join(cells) + "|"


def _format_grid(
    plot_ids: list[str],
    labels: list[str],
    n_cols: int,
    cell_width: int,
    sep: str,
) -> list[str]:
    """Render a flat list of (plot_id, label) pairs into a grid of rows.

    Returns a list of lines (no trailing newline).
    """
    assert len(plot_ids) == len(labels)
    lines: list[str] = []
    hrule = _hrule(n_cols, cell_width, sep="-")

    for chunk_start in range(0, len(plot_ids), n_cols):
        chunk_ids = plot_ids[chunk_start: chunk_start + n_cols]
        chunk_labels = labels[chunk_start: chunk_start + n_cols]

        # Pad last row to full width for alignment.
        pad = n_cols - len(chunk_ids)
        chunk_ids = chunk_ids + [""] * pad
        chunk_labels = chunk_labels + [""] * pad

        lines.append(hrule)
        lines.append(_row_line(chunk_ids, cell_width))
        lines.append(_row_line(chunk_labels, cell_width))

    lines.append(hrule)
    return lines


# ---------------------------------------------------------------------------
# CRD renderer
# ---------------------------------------------------------------------------

def _render_crd(
    design: object,
    mode: str,
    row_length: int | None,
    separator: str,
    cell_width: int,
) -> str:
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]

    header = (
        f"{meta.design_type} — {meta.n_treatments} treatments × "
        f"{meta.replications} replications = {meta.total_plots} plots"
    )
    seed_note = f"Seed: {meta.seed}" if meta.seed is not None else "Seed: None (random)"
    lines = [header, seed_note, ""]

    plot_ids: list[str] = tbl["plot_id"].tolist()
    labels: list[str] = tbl["treatment"].tolist()

    # Determine cell width: max of requested minimum and longest label/id.
    actual_width = max(
        cell_width,
        max(len(s) for s in plot_ids + labels),
    )

    if mode == "sequence":
        # One column per plot, wrapped at a reasonable terminal width.
        n_cols = row_length or min(meta.total_plots, 10)
        lines.extend(_format_grid(plot_ids, labels, n_cols, actual_width, separator))
    elif mode in ("grid", "auto"):
        n_cols = row_length or meta.extra.get("row_length") or meta.n_treatments
        lines.extend(_format_grid(plot_ids, labels, n_cols, actual_width, separator))
    else:
        raise ValueError(f"Unknown CRD mode {mode!r}. Choose 'grid', 'sequence', or 'auto'.")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# RCBD renderer
# ---------------------------------------------------------------------------

def _render_rcbd(
    design: object,
    mode: str,
    row_length: int | None,
    block_separator: bool,
    separator: str,
    cell_width: int,
) -> str:
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    n_blocks = meta.replications

    if mode not in ("blocks", "grid", "sequence"):
        raise ValueError(f"Unknown RCBD mode {mode!r}. Choose 'blocks', 'grid', or 'sequence'.")

    header = (
        f"{meta.design_type} — {meta.n_treatments} treatments × "
        f"{n_blocks} blocks = {meta.total_plots} plots"
    )
    seed_note = f"Seed: {meta.seed}" if meta.seed is not None else "Seed: None (random)"
    lines = [header, seed_note, ""]

    all_labels = tbl["treatment"].tolist() + tbl["plot_id"].tolist()
    actual_width = max(cell_width, max(len(s) for s in all_labels))
    n_cols = row_length or meta.n_treatments

    # Build a compact separator
    total_w = n_cols * (actual_width + 3) + 1
    eq_line = (separator or "=") * total_w

    blocks = sorted(tbl["block"].unique())
    for b_idx, block_id in enumerate(blocks):
        block_df = tbl[tbl["block"] == block_id].sort_values("position_in_block")
        plot_ids: list[str] = block_df["plot_id"].tolist()
        labels: list[str] = block_df["treatment"].tolist()

        lines.append(f"BLOCK {block_id}")
        lines.extend(_format_grid(plot_ids, labels, n_cols, actual_width, "-"))
        lines.append("")  # blank line after each block

        if block_separator and b_idx < len(blocks) - 1:
            lines.append(eq_line)
            lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Latin square renderer
# ---------------------------------------------------------------------------

def _render_latin_square(design: object, cell_width: int) -> str:
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    t = meta.n_treatments

    header = f"Latin square — {t} treatments, {t} × {t} = {meta.total_plots} plots"
    seed_note = f"Seed: {meta.seed}" if meta.seed is not None else "Seed: None (random)"
    lines = [header, seed_note, ""]

    width = max(cell_width, max(len(s) for s in tbl["plot_id"].tolist() + tbl["treatment"].tolist()))
    hrule = _hrule(t + 1, width, "-")
    lines.append(hrule)
    lines.append(_row_line([""] + [f"Col {c}" for c in range(1, t + 1)], width))
    for r in range(1, t + 1):
        rdf = tbl[tbl["row"] == r].sort_values("column")
        lines.append(hrule)
        lines.append(_row_line([f"Row {r}"] + rdf["plot_id"].tolist(), width))
        lines.append(_row_line([""] + rdf["treatment"].tolist(), width))
    lines.append(hrule)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Split-plot renderer
# ---------------------------------------------------------------------------

def _render_split_plot(
    design: object,
    block_separator: bool,
    separator: str,
    cell_width: int,
) -> str:
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    main_name = meta.extra["main_name"]
    sub_name = meta.extra["sub_name"]
    n_main = len(meta.extra["main_levels"])
    n_sub = len(meta.extra["sub_levels"])

    header = (
        f"Split-plot — main factor '{main_name}' ({n_main} levels) × "
        f"sub factor '{sub_name}' ({n_sub} levels), "
        f"{meta.replications} blocks = {meta.total_plots} plots"
    )
    seed_note = f"Seed: {meta.seed}" if meta.seed is not None else "Seed: None (random)"
    lines = [header, seed_note, ""]

    width = max(cell_width, max(len(s) for s in tbl["plot_id"].tolist() + tbl[sub_name].tolist()))
    total_w = n_sub * (width + 3) + 1
    blocks = sorted(tbl["block"].unique())
    for b_idx, block_id in enumerate(blocks):
        bdf = tbl[tbl["block"] == block_id]
        lines.append(f"BLOCK {block_id}")
        for mp in bdf["main_plot"].drop_duplicates().tolist():
            mdf = bdf[bdf["main_plot"] == mp]
            lines.append(f"  Main plot {mp}: {main_name} = {mdf[main_name].iloc[0]}")
            lines.extend(
                _format_grid(mdf["plot_id"].tolist(), mdf[sub_name].tolist(), n_sub, width, "-")
            )
        lines.append("")
        if block_separator and b_idx < len(blocks) - 1:
            lines.append((separator or "=") * total_w)
            lines.append("")
    return "\n".join(lines)

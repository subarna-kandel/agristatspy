"""
Split-plot design generator (main plots in an RCBD).
"""
from __future__ import annotations

import random

import pandas as pd

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.validation import (
    validate_plot_numbering,
    validate_replications,
    validate_seed,
)
from agristatspy.exceptions import DesignError


def _levels(name: str, levels: object) -> list[str]:
    if not isinstance(levels, (list, tuple)) or len(levels) < 2:
        raise DesignError(f"'{name}' needs a list of at least 2 level labels.")
    lv = [str(x) for x in levels]
    if any(not x.strip() for x in lv) or len(set(lv)) != len(lv):
        raise DesignError(f"'{name}' levels must be non-empty and unique.")
    return lv


def split_plot(
    main_levels: list[str],
    sub_levels: list[str],
    replications: int = 1,
    main_name: str = "main",
    sub_name: str = "sub",
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate a split-plot design with main plots arranged in an RCBD.

    Within each block, the levels of the main-plot factor are randomized to
    whole (main) plots.  Each main plot is then split and the levels of the
    sub-plot factor are randomized independently within it.

    Parameters
    ----------
    main_levels:
        Levels of the main-plot (whole-plot) factor, at least 2.
    sub_levels:
        Levels of the sub-plot factor, at least 2.
    replications:
        Number of complete blocks.
    main_name, sub_name:
        Column names for the two factors in the output table.
    seed, plot_prefix, start_plot:
        As for :func:`rcbd`.

    Returns
    -------
    DesignResult
        Table columns: ``plot_id``, ``block``, ``main_plot``, ``<main_name>``,
        ``<sub_name>``, ``treatment`` (combined label ``"main | sub"``).

    Raises
    ------
    DesignError
        If the arguments are invalid.
    """
    main = _levels("main_levels", main_levels)
    sub = _levels("sub_levels", sub_levels)
    n_blocks = validate_replications(replications)
    sd = validate_seed(seed)
    plot_prefix, start_plot = validate_plot_numbering(plot_prefix, start_plot)
    if (
        not isinstance(main_name, str)
        or not isinstance(sub_name, str)
        or not main_name
        or not sub_name
        or main_name == sub_name
        or {main_name, sub_name} & {"plot_id", "block", "main_plot", "treatment"}
    ):
        raise DesignError(
            "'main_name' and 'sub_name' must be distinct non-empty strings that "
            "do not clash with 'plot_id', 'block', 'main_plot', or 'treatment'."
        )

    total = len(main) * len(sub) * n_blocks
    width = max(2, len(str(start_plot + total - 1)))
    rng = random.Random(sd)

    records = []
    plot_num = start_plot
    for b in range(1, n_blocks + 1):
        order = list(main)
        rng.shuffle(order)
        for mp_idx, m in enumerate(order, start=1):
            sub_order = list(sub)
            rng.shuffle(sub_order)
            for s in sub_order:
                records.append(
                    {
                        "plot_id": f"{plot_prefix}{str(plot_num).zfill(width)}",
                        "block": b,
                        "main_plot": f"B{b}-M{mp_idx}",
                        main_name: m,
                        sub_name: s,
                        "treatment": f"{m} | {s}",
                    }
                )
                plot_num += 1

    labels = tuple(f"{m} | {s}" for m in main for s in sub)
    meta = DesignMetadata(
        design_type="SplitPlot",
        treatments=labels,
        n_treatments=len(labels),
        replications=n_blocks,
        total_plots=total,
        seed=sd,
        plot_prefix=plot_prefix,
        start_plot=start_plot,
        extra={
            "main_name": main_name,
            "sub_name": sub_name,
            "main_levels": list(main),
            "sub_levels": list(sub),
        },
    )
    result = DesignResult(table=pd.DataFrame(records), metadata=meta)
    result.validate()
    return result

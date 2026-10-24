"""
Latin Square Design generator.
"""
from __future__ import annotations

import random

import pandas as pd

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.validation import (
    parse_treatments,
    validate_plot_numbering,
    validate_seed,
)
from agristatspy.exceptions import DesignError


def latin_square(
    treatments: list[str] | None = None,
    n_treatments: int | None = None,
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate a randomized Latin Square Design.

    A Latin square of ``t`` treatments has ``t`` rows and ``t`` columns, and
    each treatment appears exactly once in every row and every column.  Rows
    and columns act as two crossed blocking factors (for example, a fertility
    gradient in one direction and moisture in the other).

    The square is built from the cyclic square and randomized by permuting
    rows, columns, and treatment labels independently.  This yields a valid,
    randomized square but does not sample uniformly from all Latin squares of
    that order.

    Parameters
    ----------
    treatments, n_treatments:
        As for :func:`agristatspy.design.crd`.  At least 2 treatments are
        required.
    seed:
        Integer seed for reproducible randomization, or ``None``.
    plot_prefix, start_plot:
        Plot-ID prefix and first plot number.

    Returns
    -------
    DesignResult
        Table columns: ``plot_id``, ``row``, ``column``, ``treatment``.
        ``metadata.replications`` equals the number of treatments.

    Raises
    ------
    DesignError
        If parameters are invalid.
    """
    labels = parse_treatments(treatments, n_treatments)
    sd = validate_seed(seed)
    plot_prefix, start_plot = validate_plot_numbering(plot_prefix, start_plot)

    t = len(labels)
    if t < 2:
        raise DesignError("A Latin square requires at least 2 treatments.")

    rng = random.Random(sd)
    row_perm = list(range(t))
    col_perm = list(range(t))
    trt_perm = list(labels)
    rng.shuffle(row_perm)
    rng.shuffle(col_perm)
    rng.shuffle(trt_perm)

    total = t * t
    width = max(2, len(str(start_plot + total - 1)))
    records = []
    plot_num = start_plot
    for r in range(t):
        for c in range(t):
            symbol = (row_perm[r] + col_perm[c]) % t
            records.append(
                {
                    "plot_id": f"{plot_prefix}{str(plot_num).zfill(width)}",
                    "row": r + 1,
                    "column": c + 1,
                    "treatment": trt_perm[symbol],
                }
            )
            plot_num += 1

    meta = DesignMetadata(
        design_type="LSQ",
        treatments=tuple(labels),
        n_treatments=t,
        replications=t,
        total_plots=total,
        seed=sd,
        plot_prefix=plot_prefix,
        start_plot=start_plot,
        extra={"size": f"{t} × {t}"},
    )
    result = DesignResult(table=pd.DataFrame(records), metadata=meta)
    result.validate()
    return result

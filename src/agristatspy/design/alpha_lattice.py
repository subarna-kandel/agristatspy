"""
Alpha-lattice (resolvable incomplete block) design generator.
"""
from __future__ import annotations

import itertools
import random

import pandas as pd

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.validation import (
    parse_treatments,
    validate_plot_numbering,
    validate_replications,
    validate_seed,
)
from agristatspy.exceptions import DesignError


def alpha_lattice(
    treatments: list[str] | None = None,
    n_treatments: int | None = None,
    block_size: int = 4,
    replications: int = 2,
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate an alpha-lattice design (resolvable incomplete blocks).

    Each replicate contains every treatment exactly once, split into
    ``s = t / block_size`` incomplete blocks of ``block_size`` plots.  The
    layout is built from the cyclic alpha(0,1) generating array of Patterson
    and Williams (the treatment in row *j* of block *c* of replicate *i* is
    ``(j * i + c) mod s + s * j``), then randomized: treatment labels,
    block order within each replicate, and plot order within each block.

    Parameters
    ----------
    treatments, n_treatments:
        Treatment labels or their number, as for :func:`rcbd`.  The number of
        treatments must be a multiple of *block_size*.
    block_size:
        Plots per incomplete block (``k >= 2``).
    replications:
        Number of replicates (``r``).  Must not exceed the number of blocks
        per replicate (``s = t / k``).
    seed, plot_prefix, start_plot:
        As for :func:`rcbd`.

    Returns
    -------
    DesignResult
        Table columns: ``plot_id``, ``replicate``, ``block`` (within
        replicate), ``position_in_block``, ``treatment``.  Metadata ``extra``
        holds ``block_size``, ``blocks_per_rep`` and the minimum and maximum
        number of times a pair of treatments shares a block.

    Raises
    ------
    DesignError
        If *treatments* is not a multiple of *block_size* or *replications*
        exceeds the blocks per replicate.

    Notes
    -----
    The cyclic generating array gives well-connected designs when ``s`` is
    prime; for composite ``s`` check ``extra['max_concurrence']`` and
    consider fewer replicates.  Analyse with :func:`agristatspy.anova_alpha_lattice`.
    """
    labels = parse_treatments(treatments, n_treatments)
    r = validate_replications(replications)
    sd = validate_seed(seed)
    plot_prefix, start_plot = validate_plot_numbering(plot_prefix, start_plot)
    k = block_size
    if isinstance(k, bool) or not isinstance(k, int) or k < 2:
        raise DesignError(f"block_size must be an integer >= 2; got {block_size!r}.")
    t = len(labels)
    if t % k:
        raise DesignError(
            f"Number of treatments ({t}) must be a multiple of block_size ({k})."
        )
    s = t // k
    if s < 2:
        raise DesignError("An alpha lattice needs at least 2 blocks per replicate.")
    if r > s:
        raise DesignError(
            f"replications ({r}) cannot exceed blocks per replicate ({s}) for the "
            "cyclic alpha(0,1) array; use fewer replications or a smaller block_size."
        )

    rng = random.Random(sd)
    names = list(labels)
    rng.shuffle(names)  # random assignment of labels to array codes 0..t-1

    records = []
    width = max(2, len(str(start_plot + t * r - 1)))
    plot_num = start_plot
    for i in range(r):
        blocks = []
        for c in range(s):
            blocks.append([names[(j * i + c) % s + s * j] for j in range(k)])
        rng.shuffle(blocks)
        for b, members in enumerate(blocks, start=1):
            rng.shuffle(members)
            for pos, trt in enumerate(members, start=1):
                records.append(
                    {
                        "plot_id": f"{plot_prefix}{str(plot_num).zfill(width)}",
                        "replicate": i + 1,
                        "block": b,
                        "position_in_block": pos,
                        "treatment": trt,
                    }
                )
                plot_num += 1
    table = pd.DataFrame(records)

    conc = dict.fromkeys(itertools.combinations(sorted(labels), 2), 0)
    for _, g in table.groupby(["replicate", "block"]):
        for pair in itertools.combinations(sorted(g["treatment"]), 2):
            conc[pair] += 1
    meta = DesignMetadata(
        design_type="AlphaLattice",
        treatments=tuple(labels),
        n_treatments=t,
        replications=r,
        total_plots=t * r,
        seed=sd,
        plot_prefix=plot_prefix,
        start_plot=start_plot,
        extra={
            "block_size": k,
            "blocks_per_rep": s,
            "min_concurrence": min(conc.values()),
            "max_concurrence": max(conc.values()),
        },
    )
    result = DesignResult(table=table, metadata=meta)
    result.validate()
    return result

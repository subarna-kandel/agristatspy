"""
Randomized Complete Block Design (RCBD/RBD) generator.
"""
from __future__ import annotations

import random

import pandas as pd

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.validation import (
    parse_treatments,
    validate_replications,
    validate_plot_numbering,
    validate_seed,
)


def rcbd(
    treatments: list[str] | None = None,
    n_treatments: int | None = None,
    replications: int = 1,
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate a Randomized Complete Block Design (RCBD/RBD).

    In a standard RCBD, each block contains every treatment exactly once.
    Treatments are independently randomized within each block.

    Exactly one of *treatments* or *n_treatments* must be supplied (or both,
    provided they are consistent).

    Parameters
    ----------
    treatments:
        List of unique treatment label strings.  If omitted, labels are
        auto-generated as ``T1``, ``T2``, … from *n_treatments*.
    n_treatments:
        Number of treatments.  Required when *treatments* is not supplied.
        Must equal ``len(treatments)`` when both are given.
    replications:
        Number of complete blocks.  In RCBD terminology each replication is a
        complete block.  Must be a positive integer (boolean values are
        rejected).
    seed:
        Integer seed for the random-number generator.  Treatments are
        shuffled independently within each block, with each block receiving
        a different derived seed to avoid identical within-block orderings.
    plot_prefix:
        String prepended to each plot number.  Default is ``"P"``.
    start_plot:
        Integer number of the first plot.  Default is ``1``.

    Returns
    -------
    DesignResult
        Object containing the assignment table, metadata, validation, text
        layout, and export methods.

    Raises
    ------
    DesignError
        If parameters are invalid or contradictory.

    Notes
    -----
    The standard additive RCBD model (treatment + block) does not separately
    estimate a treatment-by-block interaction.  If interaction is suspected,
    consider a more complex design.

    Examples
    --------
    >>> import agristatspy as ag
    >>> design = ag.design.rcbd(
    ...     treatments=["Control", "N50", "N100", "N150"],
    ...     replications=3,
    ...     seed=42,
    ... )
    >>> print(design.text_layout())
    """
    # ---- Validation ----
    labels = parse_treatments(treatments, n_treatments)
    n_blocks = validate_replications(replications)
    sd = validate_seed(seed)
    plot_prefix, start_plot = validate_plot_numbering(plot_prefix, start_plot)

    k = len(labels)
    total = k * n_blocks

    # ---- Randomization ----
    # Use a master RNG seeded from 'sd' to derive per-block seeds.
    # This ensures that block 1 does not always have the same treatment order
    # regardless of the number of blocks requested.
    master_rng = random.Random(sd)

    records = []
    plot_num = start_plot
    width = max(2, len(str(start_plot + total - 1)))

    for block_idx in range(1, n_blocks + 1):
        # Each block gets an independent shuffle.
        block_labels = list(labels)
        master_rng.shuffle(block_labels)

        for pos_idx, treatment in enumerate(block_labels, start=1):
            plot_id = f"{plot_prefix}{str(plot_num).zfill(width)}"
            records.append(
                {
                    "plot_id": plot_id,
                    "block": block_idx,
                    "position_in_block": pos_idx,
                    "treatment": treatment,
                }
            )
            plot_num += 1

    table = pd.DataFrame(records)

    meta = DesignMetadata(
        design_type="RCBD",
        treatments=tuple(labels),
        n_treatments=k,
        replications=n_blocks,
        total_plots=total,
        seed=sd,
        plot_prefix=plot_prefix,
        start_plot=start_plot,
        extra={"n_blocks": n_blocks},
    )

    result = DesignResult(table=table, metadata=meta)
    result.validate()  # Always validate on creation.
    return result

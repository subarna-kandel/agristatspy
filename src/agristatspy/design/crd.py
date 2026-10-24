"""
Completely Randomized Design (CRD) generator.
"""
from __future__ import annotations

import random

import pandas as pd

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.validation import (
    parse_treatments,
    validate_replications,
    validate_seed,
)


def crd(
    treatments: list[str] | None = None,
    n_treatments: int | None = None,
    replications: int = 1,
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate a Completely Randomized Design (CRD).

    In a CRD, every experimental unit is independently assigned to a
    treatment using simple random allocation with no blocking factor.

    Exactly one of *treatments* or *n_treatments* must be supplied (or both,
    provided they are consistent with each other).

    Parameters
    ----------
    treatments:
        List of unique treatment label strings.  If omitted, labels are
        auto-generated as ``T1``, ``T2``, … from *n_treatments*.
    n_treatments:
        Number of treatments.  Required when *treatments* is not supplied.
        Must equal ``len(treatments)`` when both are given.
    replications:
        Number of experimental units per treatment.  Must be a positive
        integer (boolean values are rejected).
    seed:
        Integer seed for the random-number generator.  Provides reproducible
        assignments within the same package version and Python environment.
        Pass ``None`` for a non-reproducible random allocation.
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

    Examples
    --------
    >>> import agristatspy as ag
    >>> design = ag.design.crd(
    ...     treatments=["Control", "N50", "N100"],
    ...     replications=4,
    ...     seed=12,
    ... )
    >>> print(design.summary())
    """
    # ---- Validation ----
    labels = parse_treatments(treatments, n_treatments)
    reps = validate_replications(replications)
    sd = validate_seed(seed)

    k = len(labels)
    total = k * reps

    # ---- Randomization ----
    rng = random.Random(sd)
    # Build a list of all (treatment, within-treatment-replicate-index) pairs
    # and shuffle them to get the randomized field order.
    pool: list[str] = labels * reps
    rng.shuffle(pool)

    # Track within-treatment replicate counters for the replication_number column.
    rep_counters: dict[str, int] = dict.fromkeys(labels, 0)
    records = []
    for i, treatment in enumerate(pool):
        rep_counters[treatment] += 1
        plot_num = start_plot + i
        width = max(2, len(str(start_plot + total - 1)))
        plot_id = f"{plot_prefix}{str(plot_num).zfill(width)}"
        records.append(
            {
                "plot_id": plot_id,
                "treatment": treatment,
                "replication_number": rep_counters[treatment],
            }
        )

    table = pd.DataFrame(records)

    meta = DesignMetadata(
        design_type="CRD",
        treatments=tuple(labels),
        n_treatments=k,
        replications=reps,
        total_plots=total,
        seed=sd,
        plot_prefix=plot_prefix,
        start_plot=start_plot,
    )

    result = DesignResult(table=table, metadata=meta)
    result.validate()  # Always validate on creation.
    return result

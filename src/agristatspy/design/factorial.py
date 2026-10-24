"""
Factorial design generator (CRD or RCBD layout).
"""
from __future__ import annotations

import itertools

from agristatspy.design.base import DesignMetadata, DesignResult
from agristatspy.design.crd import crd
from agristatspy.design.rcbd import rcbd
from agristatspy.exceptions import DesignError

_SEP = " | "


def _validate_factors(factors: object) -> dict[str, list[str]]:
    if not isinstance(factors, dict) or len(factors) < 2:
        raise DesignError(
            "'factors' must be a dict of at least two factors, "
            "e.g. {'N': ['N0', 'N50'], 'Var': ['V1', 'V2']}."
        )
    clean: dict[str, list[str]] = {}
    for name, levels in factors.items():
        if not isinstance(name, str) or not name.strip():
            raise DesignError(f"Factor names must be non-empty strings; got {name!r}.")
        if not isinstance(levels, (list, tuple)) or len(levels) < 2:
            raise DesignError(f"Factor '{name}' needs a list of at least 2 levels.")
        lv = [str(x) for x in levels]
        if any(not x.strip() for x in lv):
            raise DesignError(f"Factor '{name}' has an empty level label.")
        if len(set(lv)) != len(lv):
            raise DesignError(f"Factor '{name}' has duplicate levels.")
        if any(_SEP in x for x in lv):
            raise DesignError(f"Levels of '{name}' must not contain {_SEP!r}.")
        clean[name] = lv
    return clean


def factorial(
    factors: dict[str, list[str]],
    replications: int = 1,
    layout: str = "crd",
    seed: int | None = None,
    plot_prefix: str = "P",
    start_plot: int = 1,
) -> DesignResult:
    """Generate a factorial experiment laid out as a CRD or an RCBD.

    Every combination of factor levels is a treatment.  Combinations are
    randomized exactly as in :func:`crd` (``layout="crd"``) or :func:`rcbd`
    (``layout="rcbd"``, where *replications* is the number of blocks).

    Parameters
    ----------
    factors:
        Mapping of factor name to a list of level labels, with at least two
        factors and two levels each.
    replications:
        Replicates per combination (CRD) or number of blocks (RCBD).
    layout:
        ``"crd"`` or ``"rcbd"``.
    seed, plot_prefix, start_plot:
        As for :func:`crd`.

    Returns
    -------
    DesignResult
        The base table plus one column per factor.  The ``treatment`` column
        holds the combined label, e.g. ``"N50 | V2"``.

    Raises
    ------
    DesignError
        If the factors or layout are invalid.
    """
    fac = _validate_factors(factors)
    if layout not in ("crd", "rcbd"):
        raise DesignError(f"'layout' must be 'crd' or 'rcbd'; got {layout!r}.")

    combos = list(itertools.product(*fac.values()))
    labels = [_SEP.join(c) for c in combos]
    lookup = dict(zip(labels, combos, strict=True))

    kwargs = {
        "treatments": labels,
        "replications": replications,
        "seed": seed,
        "plot_prefix": plot_prefix,
        "start_plot": start_plot,
    }
    base = crd(**kwargs) if layout == "crd" else rcbd(**kwargs)  # type: ignore[arg-type]

    table = base.table
    for i, name in enumerate(fac):
        table[name] = table["treatment"].map(lambda t, i=i: lookup[t][i])

    bm = base.metadata
    extra = dict(bm.extra)
    extra["base_design"] = bm.design_type
    extra["factors"] = {k: list(v) for k, v in fac.items()}
    meta = DesignMetadata(
        design_type=f"Factorial-{bm.design_type}",
        treatments=bm.treatments,
        n_treatments=bm.n_treatments,
        replications=bm.replications,
        total_plots=bm.total_plots,
        seed=bm.seed,
        plot_prefix=bm.plot_prefix,
        start_plot=bm.start_plot,
        extra=extra,
    )
    result = DesignResult(table=table, metadata=meta)
    result.validate()
    return result

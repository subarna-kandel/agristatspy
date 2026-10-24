"""
Design-parameter validation utilities shared by CRD and RCBD generators.
"""
from __future__ import annotations

from agristatspy.exceptions import DesignError


def parse_treatments(
    treatments: list[str] | None,
    n_treatments: int | None,
) -> list[str]:
    """Validate and resolve the treatment-label list.

    Exactly one of *treatments* or *n_treatments* must be supplied, or both
    provided they are consistent.  Returns a validated list of unique label
    strings.

    Parameters
    ----------
    treatments:
        Explicit list of treatment labels, or ``None``.
    n_treatments:
        Number of treatments when labels are not supplied, or ``None``.

    Returns
    -------
    list[str]
        Validated list of treatment-label strings.

    Raises
    ------
    DesignError
        On any invalid or contradictory argument.
    """
    if treatments is None and n_treatments is None:
        raise DesignError(
            "Supply either 'treatments' (a list of labels) or "
            "'n_treatments' (a positive integer)."
        )

    if n_treatments is not None:
        # Reject bool masquerading as int: bool is a subclass of int in Python.
        if isinstance(n_treatments, bool):
            raise DesignError("'n_treatments' must be an integer, not bool.")
        if not isinstance(n_treatments, int) or n_treatments < 1:
            raise DesignError(
                f"'n_treatments' must be a positive integer; got {n_treatments!r}."
            )

    if treatments is not None:
        if not isinstance(treatments, list) or len(treatments) == 0:
            raise DesignError(
                "'treatments' must be a non-empty list of label strings."
            )
        for label in treatments:
            if not isinstance(label, str) or label.strip() == "":
                raise DesignError(
                    f"Each treatment label must be a non-empty string; got {label!r}."
                )
        if len(treatments) != len(set(treatments)):
            duplicates = [t for t in treatments if treatments.count(t) > 1]
            raise DesignError(
                f"Treatment labels must be unique; duplicates found: {sorted(set(duplicates))}."
            )

        if n_treatments is not None and n_treatments != len(treatments):
            raise DesignError(
                f"'n_treatments' ({n_treatments}) does not match the number of "
                f"unique labels supplied ({len(treatments)})."
            )

        return list(treatments)

    # Auto-generate labels: T1, T2, …
    return [f"T{i}" for i in range(1, n_treatments + 1)]  # type: ignore[arg-type]


def validate_replications(replications: object) -> int:
    """Validate the replications argument.

    Parameters
    ----------
    replications:
        Must be a positive (non-bool) integer.

    Returns
    -------
    int
        The validated replication count.

    Raises
    ------
    DesignError
        If the value is invalid.
    """
    if isinstance(replications, bool):
        raise DesignError("'replications' must be an integer, not bool.")
    if not isinstance(replications, int) or replications < 1:
        raise DesignError(
            f"'replications' must be a positive integer; got {replications!r}."
        )
    return replications


def validate_seed(seed: object) -> int | None:
    """Validate the random seed.

    Parameters
    ----------
    seed:
        Must be an integer or ``None``.

    Returns
    -------
    int or None

    Raises
    ------
    DesignError
        If the value is not a valid seed.
    """
    if seed is None:
        return None
    if isinstance(seed, bool):
        raise DesignError("'seed' must be an integer, not bool.")
    if not isinstance(seed, int):
        raise DesignError(
            f"'seed' must be an integer or None; got {type(seed).__name__!r}."
        )
    return seed

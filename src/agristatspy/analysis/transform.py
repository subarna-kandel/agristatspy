"""
Data transformations for meeting ANOVA assumptions.

A transformation is applied to the response *before* the analysis and the
original column is never modified or dropped.  Missing values stay missing.
Invalid input (for example ``log`` of a non-positive value) raises an
:class:`~agristatspy.exceptions.AnalysisError` instead of being silently
dropped or shifted.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.exceptions import AnalysisError

METHODS = (
    "log",
    "log10",
    "sqrt",
    "arcsine",
    "reciprocal",
    "square",
    "boxcox",
)


@dataclass
class TransformResult:
    """A transformed copy of the data and what is needed to undo it."""

    data: pd.DataFrame
    column: str
    source: str
    method: str
    params: dict = field(default_factory=dict)
    note: str = ""

    def back_transform(self, values):
        """Return ``values`` (e.g. treatment means) on the original scale.

        Back-transformed means are not arithmetic means of the original data
        (for ``log`` they estimate the geometric mean); report them as such.
        """
        y = np.asarray(values, dtype=float)
        m, p = self.method, self.params
        if m == "log":
            out = np.exp(y) - p["constant"]
        elif m == "log10":
            out = 10.0**y - p["constant"]
        elif m == "sqrt":
            out = y**2 - p["constant"]
        elif m == "arcsine":
            out = np.sin(np.radians(y)) ** 2
            if p["scale"] == "percent":
                out = out * 100.0
        elif m == "reciprocal":
            out = 1.0 / y - p["constant"]
        elif m == "square":
            out = np.sqrt(y)
        else:  # boxcox
            lam = p["lambda"]
            out = np.exp(y) if lam == 0 else np.power(lam * y + 1.0, 1.0 / lam)
        return out if np.ndim(values) else float(out)

    def __str__(self) -> str:
        pars = ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                         for k, v in self.params.items())
        head = f"Transformation: {self.method}({self.source}) -> {self.column}"
        return f"{head}\n  parameters: {pars or 'none'}\n  {self.note}"


def _note(method: str) -> str:
    return {
        "log": "Back-transformed means estimate geometric means.",
        "log10": "Back-transformed means estimate geometric means.",
        "sqrt": "Typical for counts (use constant=0.5 when zeros occur).",
        "arcsine": "For proportions/percentages; not suitable for values "
                   "mostly within 30-70% (little gain) or exactly 0/100 "
                   "without care.",
        "reciprocal": "Reverses the order of values; interpret with care.",
        "square": "For left-skewed data.",
        "boxcox": "Lambda estimated by maximum likelihood on the pooled "
                  "response; re-check residuals after transforming.",
    }[method]


def transform(
    data: pd.DataFrame,
    column: str,
    method: str,
    constant: float = 0.0,
    scale: str = "percent",
    lmbda: float | None = None,
    name: str | None = None,
) -> TransformResult:
    """Transform a response column and return a copy of the data.

    Parameters
    ----------
    data:
        Tidy DataFrame.  It is not modified.
    column:
        Numeric column to transform.
    method:
        ``"log"`` (natural log of ``x + constant``), ``"log10"``,
        ``"sqrt"`` (``sqrt(x + constant)``; use ``constant=0.5`` for counts
        with zeros), ``"arcsine"`` (``arcsin(sqrt(p))`` in degrees, for
        proportions or percentages), ``"reciprocal"`` (``1 / (x + constant)``),
        ``"square"`` or ``"boxcox"`` (requires ``x > 0``).
    constant:
        Value added before ``log``/``log10``/``sqrt``/``reciprocal``.  It is
        never chosen automatically.
    scale:
        For ``"arcsine"`` only: ``"percent"`` (0-100) or ``"proportion"`` (0-1).
    lmbda:
        For ``"boxcox"`` only: use this lambda instead of estimating it.
    name:
        Name of the new column (default ``"<column>_<method>"``).

    Raises
    ------
    AnalysisError
        On an unknown method or column, a non-numeric column, values outside
        the domain of the transformation, or an existing target column.
    """
    if method not in METHODS:
        raise AnalysisError(f"Unknown transformation {method!r}; choose one of {METHODS}.")
    if column not in data.columns:
        raise AnalysisError(f"Column {column!r} not found in data.")
    if not pd.api.types.is_numeric_dtype(data[column]):
        raise AnalysisError(f"Column {column!r} must be numeric.")
    if scale not in ("percent", "proportion"):
        raise AnalysisError("scale must be 'percent' or 'proportion'.")
    new = name or f"{column}_{method}"
    if new in data.columns:
        raise AnalysisError(f"Column {new!r} already exists; pass name= to choose another.")

    x = data[column].astype(float)
    obs = x.dropna()
    if obs.empty:
        raise AnalysisError(f"Column {column!r} has no observed values.")
    if not np.isfinite(obs).all():
        raise AnalysisError(f"Column {column!r} contains infinite values.")

    params: dict = {}
    if method in ("log", "log10", "sqrt", "reciprocal"):
        shifted = obs + constant
        params["constant"] = float(constant)
        if method == "sqrt":
            if (shifted < 0).any():
                raise AnalysisError("sqrt requires x + constant >= 0.")
        elif (shifted <= 0).any():
            raise AnalysisError(
                f"{method} requires x + constant > 0 (found {int((shifted <= 0).sum())} "
                "non-positive value(s)); set constant= explicitly (e.g. 1 for counts "
                "with zeros)."
            )
        y = {
            "log": np.log,
            "log10": np.log10,
            "sqrt": np.sqrt,
            "reciprocal": lambda v: 1.0 / v,
        }[method](x + constant)
    elif method == "arcsine":
        top = 100.0 if scale == "percent" else 1.0
        if (obs < 0).any() or (obs > top).any():
            raise AnalysisError(f"arcsine with scale={scale!r} requires values in [0, {top:g}].")
        params["scale"] = scale
        y = np.degrees(np.arcsin(np.sqrt(x / top)))
    elif method == "square":
        y = x**2
    else:  # boxcox
        if (obs <= 0).any():
            raise AnalysisError("boxcox requires all values > 0; add a constant first.")
        if lmbda is None:
            _, lam = stats.boxcox(obs.to_numpy())
        else:
            lam = float(lmbda)
        lam = float(lam)
        params["lambda"] = lam
        y = np.log(x) if lam == 0 else (x**lam - 1.0) / lam

    out = data.copy()
    out[new] = y
    return TransformResult(
        data=out, column=new, source=column, method=method, params=params, note=_note(method)
    )

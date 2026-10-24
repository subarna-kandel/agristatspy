"""
ANOVA for Latin square, factorial, and split-plot designs, plus a dispatcher.

Models (fixed effects, balanced data)
-------------------------------------
Latin square
    y_ijk = μ + ρ_i + γ_j + τ_k + ε_ijk
Factorial (CRD)
    y = μ + all main effects and interactions of the factors + ε
Factorial (RCBD)
    as above plus a block effect (no block × treatment interaction)
Split-plot (main plots in an RCBD)
    y_ijk = μ + β_i + α_j + (βα)_ij + γ_k + (αγ)_jk + ε_ijk
    where (βα)_ij is the main-plot error (error a) and ε is the sub-plot
    error (error b).  The main-plot factor is tested against error a; the
    sub-plot factor and the interaction against error b.

Balanced layouts are required because sequential (Type I) sums of squares
are only order-independent when the data are balanced.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

from agristatspy.analysis.anova import ANOVAResult, anova_crd, anova_rcbd
from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError

_ASSUMPTIONS = (
    "Assumptions: independent experimental units, approximately normal "
    "residuals,\n  homogeneous variances.\n"
    "Statistical significance does not imply agronomic importance."
)


def _prepare(
    data: pd.DataFrame,
    response: str,
    factors: dict[str, str],
    context: str,
) -> tuple[pd.DataFrame, int]:
    """Subset, drop missing responses, and rename factors to safe codes.

    *factors* maps safe code -> user column name.
    """
    cols = [response, *factors.values()]
    _check_columns(data, cols, context)
    _check_numeric(data, response, context)
    df = data[cols].copy()
    n_before = len(df)
    df = df.dropna(subset=[response])
    missing = n_before - len(df)
    renamed = {user: code for code, user in factors.items()}
    df = df.rename(columns={**renamed, response: "y"})
    for code in factors:
        df[code] = df[code].astype(str)
    df["y"] = df["y"].astype(float)
    return df, missing


def _fit(df: pd.DataFrame, formula: str) -> pd.DataFrame:
    model = ols(formula, data=df).fit()
    return anova_lm(model, typ=1)


def _row(source: str, r: pd.Series | None, ms_err: float | None, df_err: int | None):
    """Build a table row from a statsmodels row, testing against *ms_err*."""
    from scipy import stats

    ss = float(r["sum_sq"])
    dfx = int(round(float(r["df"])))
    ms = ss / dfx
    if ms_err is None or df_err is None:
        f, p = float("nan"), float("nan")
    else:
        f = ms / ms_err
        p = float(stats.f.sf(f, dfx, df_err))
    return {"source": source, "df": dfx, "ss": ss, "ms": ms, "f_stat": f, "p_value": p}


def _error_row(source: str, ss: float, dfe: int) -> dict:
    return {
        "source": source, "df": dfe, "ss": ss, "ms": ss / dfe,
        "f_stat": float("nan"), "p_value": float("nan"),
    }


def _total_row(df: pd.DataFrame) -> dict:
    y = df["y"].values
    return {
        "source": "Total", "df": len(y) - 1,
        "ss": float(np.sum((y - y.mean()) ** 2)),
        "ms": float("nan"), "f_stat": float("nan"), "p_value": float("nan"),
    }


def _cv(ms_err: float, grand_mean: float) -> float:
    return float(np.sqrt(ms_err) / abs(grand_mean) * 100) if abs(grand_mean) > 1e-10 else float("nan")


def _require_unique_cells(df: pd.DataFrame, keys: list[str], expected: int, context: str) -> None:
    sizes = df.groupby(keys).size()
    if (sizes != 1).any() or len(sizes) != expected:
        raise AnalysisError(
            f"{context} requires exactly one observation for each combination of "
            f"{keys} (expected {expected} combinations, found {len(sizes)}); "
            "missing or duplicated combinations are not allowed."
        )


# ---------------------------------------------------------------------------
# Latin square
# ---------------------------------------------------------------------------

def anova_latin_square(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    row: str,
    column: str,
) -> ANOVAResult:
    """Additive ANOVA for a Latin Square Design (rows + columns + treatments).

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per plot.
    response, treatment, row, column:
        Column names.

    Returns
    -------
    ANOVAResult

    Raises
    ------
    AnalysisError
        If columns are missing, the response is non-numeric, the square is not
        t × t with t treatments, or any plot is missing/duplicated.

    Notes
    -----
    The model assumes no interactions among rows, columns, and treatments.
    A Latin square gives few error degrees of freedom for small t
    (``(t-1)(t-2)``), so power is low for t < 5.
    """
    ctx = "anova_latin_square"
    df, missing = _prepare(
        data, response, {"trt": treatment, "row": row, "col": column}, ctx
    )
    t = df["trt"].nunique()
    if t < 3:
        raise AnalysisError(f"{ctx} requires at least 3 treatments; found {t}.")
    if df["row"].nunique() != t or df["col"].nunique() != t:
        raise AnalysisError(
            f"{ctx} requires a square layout: {t} treatments, "
            f"{df['row'].nunique()} rows, {df['col'].nunique()} columns."
        )
    _require_unique_cells(df, ["row", "col"], t * t, ctx)
    for key in ("row", "col"):
        if (df.groupby([key, "trt"]).size() != 1).any():
            raise AnalysisError(
                f"{ctx}: each treatment must occur exactly once in every {key}."
            )

    tbl = _fit(df, "y ~ C(row) + C(col) + C(trt)")
    err = tbl.loc["Residual"]
    ms_err = float(err["sum_sq"]) / int(err["df"])
    dfe = int(err["df"])

    rows = [
        _row("Treatment", tbl.loc["C(trt)"], ms_err, dfe),
        _row("Row", tbl.loc["C(row)"], ms_err, dfe),
        _row("Column", tbl.loc["C(col)"], ms_err, dfe),
        _error_row("Error", float(err["sum_sq"]), dfe),
        _total_row(df),
    ]
    gm = float(df["y"].mean())
    notes = (
        "Model: additive fixed-effects ANOVA for a Latin square.\n"
        "  y_ijk = μ + ρ_i + γ_j + τ_k + ε_ijk\n"
        "Assumes no row × column × treatment interactions.\n"
        f"{_ASSUMPTIONS}\nCV% = sqrt(MS_error) / grand_mean × 100."
    )
    return ANOVAResult(pd.DataFrame(rows), "Latin square", len(df), gm, _cv(ms_err, gm), notes, missing)


# ---------------------------------------------------------------------------
# Factorial
# ---------------------------------------------------------------------------

def anova_factorial(
    data: pd.DataFrame,
    response: str,
    factors: list[str],
    block: str | None = None,
) -> ANOVAResult:
    """ANOVA for a balanced factorial experiment with all interactions.

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per experimental unit.
    response:
        Numeric response column.
    factors:
        Two or more factor columns.
    block:
        Optional block column for a factorial laid out as an RCBD.  Without
        it the layout is treated as a CRD.

    Returns
    -------
    ANOVAResult
        Rows for each main effect and interaction (``A``, ``B``, ``A:B``, …),
        the block (if any), Error, and Total.

    Raises
    ------
    AnalysisError
        If inputs are invalid, the data are unbalanced, or there are no error
        degrees of freedom (e.g. a single replicate without blocks).

    Notes
    -----
    With a block, the model is additive in blocks (no block × treatment
    interaction).  Higher-order interactions are tested against the pooled
    error; if they are of no interest, interpret them with care.
    """
    ctx = "anova_factorial"
    if isinstance(factors, str) or len(factors) < 2:
        raise AnalysisError(f"{ctx} requires a list of at least two factor columns.")
    if len(set(factors)) != len(factors):
        raise AnalysisError(f"{ctx}: factor columns must be distinct.")
    codes = {f"F{i}": name for i, name in enumerate(factors)}
    spec = dict(codes)
    if block is not None:
        if block in factors:
            raise AnalysisError(f"{ctx}: block column cannot also be a factor.")
        spec["blk"] = block
    df, missing = _prepare(data, response, spec, ctx)

    cell_keys = list(codes)
    sizes = df.groupby(cell_keys).size()
    n_cells = int(np.prod([df[c].nunique() for c in cell_keys]))
    if len(sizes) != n_cells or sizes.nunique() != 1:
        raise AnalysisError(
            f"{ctx} requires a balanced layout with equal replication in every "
            "factor-level combination (no missing combinations or responses)."
        )
    if block is not None:
        _require_unique_cells(df, [*cell_keys, "blk"], n_cells * df["blk"].nunique(), ctx)

    n_blocks = df["blk"].nunique() if block is not None else 1
    if len(df) - n_cells - (n_blocks - 1) < 1:
        raise AnalysisError(
            f"{ctx}: no error degrees of freedom — add replication or blocks."
        )

    rhs = " * ".join(f"C({c})" for c in codes)
    if block is not None:
        rhs = f"C(blk) + {rhs}"
    tbl = _fit(df, f"y ~ {rhs}")
    err = tbl.loc["Residual"]
    dfe = int(err["df"])
    if dfe < 1:
        raise AnalysisError(
            f"{ctx}: no error degrees of freedom — add replication or blocks."
        )
    ms_err = float(err["sum_sq"]) / dfe

    def pretty(idx: str) -> str:
        parts = []
        for term in idx.split(":"):
            code = term[2:-1]
            parts.append(codes.get(code, "Block"))
        return ":".join(parts)

    rows = []
    for idx in tbl.index:
        if idx == "Residual":
            continue
        label = "Block" if idx == "C(blk)" else pretty(idx)
        rows.append(_row(label, tbl.loc[idx], ms_err, dfe))
    rows.append(_error_row("Error", float(err["sum_sq"]), dfe))
    rows.append(_total_row(df))

    gm = float(df["y"].mean())
    layout = "RCBD" if block else "CRD"
    notes = (
        f"Model: fixed-effects factorial ANOVA ({layout} layout) with all "
        "interactions among " + ", ".join(factors) + ".\n"
        + ("Block effect is additive (no block × treatment interaction).\n" if block else "")
        + "If an interaction is significant, interpret main effects cautiously;\n"
        "  compare cell (combination) means or simple effects instead.\n"
        f"{_ASSUMPTIONS}\nCV% = sqrt(MS_error) / grand_mean × 100."
    )
    return ANOVAResult(
        pd.DataFrame(rows), f"Factorial ({layout})", len(df), gm, _cv(ms_err, gm), notes, missing
    )


# ---------------------------------------------------------------------------
# Split-plot
# ---------------------------------------------------------------------------

def anova_split_plot(
    data: pd.DataFrame,
    response: str,
    main: str,
    sub: str,
    block: str,
) -> ANOVAResult:
    """ANOVA for a split-plot design with main plots in an RCBD.

    The main-plot factor is tested against the main-plot error (error a,
    block × main interaction); the sub-plot factor and the main × sub
    interaction are tested against the sub-plot error (error b).

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per sub-plot.
    response:
        Numeric response column.
    main, sub, block:
        Main-plot factor, sub-plot factor, and block columns.

    Returns
    -------
    ANOVAResult
        ``cv_pct`` is CV(b) (from error b); CV(a) is reported in the notes.

    Raises
    ------
    AnalysisError
        If inputs are invalid or the layout is incomplete or unbalanced.
    """
    ctx = "anova_split_plot"
    if len({main, sub, block}) != 3:
        raise AnalysisError(f"{ctx}: main, sub, and block must be distinct columns.")
    df, missing = _prepare(data, response, {"A": main, "B": sub, "blk": block}, ctx)
    a, b, r = df["A"].nunique(), df["B"].nunique(), df["blk"].nunique()
    if a < 2 or b < 2 or r < 2:
        raise AnalysisError(f"{ctx} requires at least 2 main levels, 2 sub levels, and 2 blocks.")
    _require_unique_cells(df, ["blk", "A", "B"], a * b * r, ctx)

    tbl = _fit(df, "y ~ C(blk) + C(A) + C(blk):C(A) + C(B) + C(A):C(B)")
    ea = tbl.loc["C(blk):C(A)"]
    eb = tbl.loc["Residual"]
    ms_a, df_a = float(ea["sum_sq"]) / int(ea["df"]), int(ea["df"])
    ms_b, df_b = float(eb["sum_sq"]) / int(eb["df"]), int(eb["df"])

    rows = [
        _row("Block", tbl.loc["C(blk)"], None, None),
        _row(main, tbl.loc["C(A)"], ms_a, df_a),
        _error_row("Error (a)", float(ea["sum_sq"]), df_a),
        _row(sub, tbl.loc["C(B)"], ms_b, df_b),
        _row(f"{main}:{sub}", tbl.loc["C(A):C(B)"], ms_b, df_b),
        _error_row("Error (b)", float(eb["sum_sq"]), df_b),
        _total_row(df),
    ]
    gm = float(df["y"].mean())
    notes = (
        "Model: split-plot ANOVA, main plots in an RCBD.\n"
        f"  '{main}' (main plot) is tested against Error (a) = block × {main}.\n"
        f"  '{sub}' and '{main}:{sub}' are tested against Error (b).\n"
        f"{_ASSUMPTIONS}\n"
        f"CV(a)% = {_cv(ms_a, gm):.2f}  (from Error a);  "
        f"CV(b)% = {_cv(ms_b, gm):.2f}  (from Error b, reported as CV%).\n"
        "Block is not tested (treated as a blocking factor)."
    )
    return ANOVAResult(pd.DataFrame(rows), "Split-plot", len(df), gm, _cv(ms_b, gm), notes, missing)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def anova(
    data: pd.DataFrame,
    response: str,
    design: str,
    treatment: str | None = None,
    block: str | None = None,
    row: str | None = None,
    column: str | None = None,
    factors: list[str] | None = None,
    main: str | None = None,
    sub: str | None = None,
) -> ANOVAResult:
    """Run the ANOVA appropriate to *design*.

    Parameters
    ----------
    data, response:
        Tidy DataFrame and numeric response column.
    design:
        ``"crd"``, ``"rcbd"``, ``"latin_square"``, ``"factorial"``, or
        ``"split_plot"``.
    treatment, block, row, column, factors, main, sub:
        Column names; which are required depends on *design*:

        - crd: ``treatment``
        - rcbd: ``treatment``, ``block``
        - latin_square: ``treatment``, ``row``, ``column``
        - factorial: ``factors`` (optionally ``block``)
        - split_plot: ``main``, ``sub``, ``block``

    Returns
    -------
    ANOVAResult

    Raises
    ------
    AnalysisError
        If *design* is unknown or required columns are not supplied.
    """
    key = str(design).lower().replace("-", "_").replace(" ", "_")

    def need(**kw: object) -> None:
        missing = [k for k, v in kw.items() if v is None]
        if missing:
            raise AnalysisError(f"anova(design={design!r}) requires: {', '.join(missing)}.")

    if key == "crd":
        need(treatment=treatment)
        return anova_crd(data, response, treatment)  # type: ignore[arg-type]
    if key in ("rcbd", "rbd"):
        need(treatment=treatment, block=block)
        return anova_rcbd(data, response, treatment, block)  # type: ignore[arg-type]
    if key in ("latin_square", "lsq", "latin"):
        need(treatment=treatment, row=row, column=column)
        return anova_latin_square(data, response, treatment, row, column)  # type: ignore[arg-type]
    if key == "factorial":
        need(factors=factors)
        return anova_factorial(data, response, factors, block)  # type: ignore[arg-type]
    if key in ("split_plot", "splitplot"):
        need(main=main, sub=sub, block=block)
        return anova_split_plot(data, response, main, sub, block)  # type: ignore[arg-type]
    raise AnalysisError(
        f"Unknown design {design!r}. Choose from: crd, rcbd, latin_square, "
        "factorial, split_plot."
    )

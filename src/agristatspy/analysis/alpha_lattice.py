"""
Intra-block analysis of alpha-lattice (resolvable incomplete block) designs.

Model (blocks fixed, nested in replicates)
------------------------------------------
    y_ijk = μ + ρ_i + β_j(i) + τ_k + ε_ijk

Treatments are adjusted for incomplete blocks, so the treatment sum of
squares is the one *eliminating* replicates and blocks.  Adjusted means are
least-squares means; their standard errors and the standard error of every
pairwise difference come from the generalized inverse of the normal
equations, so unbalanced or incompletely connected data are handled.
Missing responses are allowed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis._export import export_tables
from agristatspy.analysis.anova import ANOVAResult
from agristatspy.analysis.comparison import ComparisonResult, compare_means
from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


def _add_letters(out: pd.DataFrame, comparisons: dict[str, ComparisonResult]) -> pd.DataFrame:
    """Append one letter-group column per comparison to a ``Treatment`` table."""
    for comp in comparisons.values():
        lt = comp.letters().set_index("treatment")["group"]
        col = "Group" if len(comparisons) == 1 else f"Group ({comp.method})"
        out[col] = out["Treatment"].map(lt).to_numpy()
    return out


def _design_matrix(codes: list[pd.Series]) -> np.ndarray:
    cols = [np.ones((len(codes[0]), 1))] if codes else []
    for c in codes:
        cols.append(pd.get_dummies(c, dtype=float).to_numpy())
    return np.hstack(cols)


def _sse(y: np.ndarray, x: np.ndarray) -> tuple[float, int]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    r = y - x @ beta
    return float(r @ r), int(np.linalg.matrix_rank(x))


@dataclass
class AlphaLatticeResult:
    """Result of :func:`anova_alpha_lattice`.

    Attributes
    ----------
    anova:
        :class:`ANOVAResult` (replicate, block within replicate, treatment
        adjusted for blocks, error, total).
    means:
        Per-treatment table: ``treatment``, ``adj_mean``, ``se``, ``raw_mean``,
        ``n``.
    sed:
        Matrix of standard errors of differences between adjusted means.
    mse, df_error:
        Intra-block error mean square and degrees of freedom.
    relative_efficiency:
        ``MSE(RCBD) / MSE(alpha lattice) × 100``: how much the incomplete-block
        analysis improves on ignoring the blocks.  Above 100 means the lattice
        helped.
    """

    anova: ANOVAResult
    means: pd.DataFrame
    sed: pd.DataFrame
    mse: float
    df_error: int
    relative_efficiency: float
    alpha: float = 0.05
    decimals: int = 2
    comparisons: dict[str, ComparisonResult] = field(default_factory=dict)
    messages: list[str] = field(default_factory=list)
    _pair_var: pd.DataFrame = field(default_factory=pd.DataFrame, repr=False)

    @property
    def avg_sed(self) -> float:
        """Average standard error of a difference between two adjusted means."""
        v = self.sed.to_numpy()
        return float(v[np.triu_indices_from(v, k=1)].mean())

    @property
    def cd(self) -> float:
        """Critical difference at ``alpha`` using the average SED."""
        return float(stats.t.ppf(1 - self.alpha / 2, self.df_error) * self.avg_sed)

    def compare(self, method: str = "tukey", alpha: float | None = None) -> ComparisonResult:
        """Pairwise comparison of adjusted means with each pair's own SED.

        *method* is ``"tukey"``, ``"lsd"``, ``"dmrt"`` or ``"bonferroni"``;
        ``.letters()`` on the result gives compact letter groups.
        """
        pv = self._pair_var
        return compare_means(
            self.means.set_index("treatment")["adj_mean"], 1, self.mse, self.df_error,
            method=method, alpha=self.alpha if alpha is None else alpha,
            pair_factor=lambda a, b: float(pv.loc[a, b]),
            error_note="intra-block error of the alpha lattice (adjusted means)",
        )

    def thesis_table(self, decimals: int | None = None) -> pd.DataFrame:
        d = self.decimals if decimals is None else decimals
        m = self.means
        out = pd.DataFrame(
            {
                "Treatment": m["treatment"].astype(str),
                "Adjusted mean ± SE": [
                    f"{a:.{d}f} ± {s:.{d}f}" for a, s in zip(m["adj_mean"], m["se"], strict=True)
                ],
                "Raw mean": m["raw_mean"].round(d),
            }
        )
        return _add_letters(out, self.comparisons)

    def tables(self) -> dict[str, pd.DataFrame]:
        t = self.anova.table.copy()
        pairwise = {f"Pairwise {c.method}"[:31]: c.table for c in self.comparisons.values()}
        return {**pairwise,
            "Adjusted means": self.thesis_table(),
            "Means (numeric)": self.means,
            "ANOVA": t,
            "Summary": pd.DataFrame(
                {
                    "statistic": ["Grand mean", "CV (%)", "Error MS", "Error df",
                                  "Average SED", f"CD ({int(self.alpha * 100)}%)",
                                  "Relative efficiency vs RCBD (%)"],
                    "value": [self.anova.grand_mean, self.anova.cv_pct, self.mse,
                              self.df_error, self.avg_sed, self.cd, self.relative_efficiency],
                }
            ),
        }

    def export(self, path: str) -> object:
        """Export to ``.xlsx``, ``.csv``, ``.md``, ``.html`` or ``.tex``."""
        return export_tables(self.tables(), path, "Alpha-lattice analysis", str(self))

    def __str__(self) -> str:
        names = ", ".join(c.method for c in self.comparisons.values())
        return "\n".join(
            [
                str(self.anova),
                "",
                f"Average SED = {self.avg_sed:.4f}   CD({int(self.alpha * 100)}%) = {self.cd:.4f}"
                f"   Relative efficiency vs RCBD = {self.relative_efficiency:.1f}%",
                "",
                "Adjusted treatment means" + (f" (letters: {names})" if names else "") + ":",
                self.thesis_table().to_string(index=False),
                *(["", *[f"• {m}" for m in self.messages]] if self.messages else []),
            ]
        )


def anova_alpha_lattice(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    replicate: str,
    block: str,
    alpha: float = 0.05,
    decimals: int = 2,
) -> AlphaLatticeResult:
    """Intra-block ANOVA and adjusted means for an alpha-lattice design.

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per plot.
    response, treatment:
        Response and treatment columns.
    replicate:
        Replicate column.
    block:
        Incomplete-block column.  Block labels may repeat across replicates
        (block 1 of replicate 1 and of replicate 2 are different blocks).
    alpha:
        Level for the critical difference.
    decimals:
        Decimals in the formatted tables.

    Returns
    -------
    AlphaLatticeResult

    Raises
    ------
    AnalysisError
        If columns are missing, the response is non-numeric, treatments are
        not connected through the blocks, or no error degrees of freedom
        remain.

    Notes
    -----
    Blocks are treated as fixed (intra-block analysis).  Recovery of
    inter-block information, as in a mixed-model REML analysis, is not done.
    Statistical significance does not establish agronomic importance.
    """
    _check_columns(data, [response, treatment, replicate, block], "anova_alpha_lattice")
    _check_numeric(data, response, "anova_alpha_lattice")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
        raise AnalysisError(f"alpha must be a number strictly between 0 and 1; got {alpha!r}.")

    df = data[[response, treatment, replicate, block]].copy()
    n_before = len(df)
    df = df.dropna(subset=[response]).reset_index(drop=True)
    missing = n_before - len(df)
    df[treatment] = df[treatment].astype(str)
    df["__rep__"] = df[replicate].astype(str)
    df["__blk__"] = df["__rep__"] + "/" + df[block].astype(str)
    y = df[response].to_numpy(float)
    n = len(y)

    x0 = np.ones((n, 1))
    x_rep = _design_matrix([df["__rep__"]])
    x_blk = _design_matrix([df["__blk__"]])  # spans replicates
    x_trt_only = _design_matrix([df["__rep__"], df[treatment]])
    x_full = np.hstack([x_blk, pd.get_dummies(df[treatment], dtype=float).to_numpy()])

    sse0, _ = _sse(y, x0)
    sse_rep, r_rep = _sse(y, x_rep)
    sse_blk, r_blk = _sse(y, x_blk)
    sse_full, r_full = _sse(y, x_full)
    sse_rcbd, r_rcbd = _sse(y, x_trt_only)

    n_trt = df[treatment].nunique()
    df_rep = r_rep - 1
    df_blk = r_blk - r_rep
    df_trt = r_full - r_blk
    df_err = n - r_full
    if df_trt < n_trt - 1:
        raise AnalysisError(
            "anova_alpha_lattice: treatments are not all connected through the blocks "
            f"(estimable treatment df = {df_trt}, expected {n_trt - 1})."
        )
    if df_err < 1:
        raise AnalysisError("anova_alpha_lattice: no error degrees of freedom; more plots needed.")

    ss_rep, ss_blk, ss_trt = sse0 - sse_rep, sse_rep - sse_blk, sse_blk - sse_full
    mse = sse_full / df_err
    ms_trt = ss_trt / df_trt
    f_trt = ms_trt / mse
    p_trt = float(stats.f.sf(f_trt, df_trt, df_err))
    nan = float("nan")
    rows = [
        ("Replicate", df_rep, ss_rep), ("Block (Replicate)", df_blk, ss_blk),
    ]
    table = pd.DataFrame(
        [{"source": s, "df": d, "ss": ss, "ms": ss / d if d else nan,
          "f_stat": nan, "p_value": nan} for s, d, ss in rows]
        + [
            {"source": "Treatment (adj.)", "df": df_trt, "ss": ss_trt, "ms": ms_trt,
             "f_stat": f_trt, "p_value": p_trt},
            {"source": "Error", "df": df_err, "ss": sse_full, "ms": mse,
             "f_stat": nan, "p_value": nan},
            {"source": "Total", "df": n - 1, "ss": sse0, "ms": nan,
             "f_stat": nan, "p_value": nan},
        ]
    )

    # Least-squares means: average design row over observations, treatment set to k.
    labels = sorted(df[treatment].unique())
    t_cols = pd.get_dummies(df[treatment], dtype=float).columns.tolist()
    nb = x_blk.shape[1]
    base = x_full[:, :nb].mean(axis=0)
    lmat = np.zeros((len(labels), x_full.shape[1]))
    for i, lab in enumerate(labels):
        lmat[i, :nb] = base
        lmat[i, nb + t_cols.index(lab)] = 1.0
    xtx_pinv = np.linalg.pinv(x_full.T @ x_full)
    beta = xtx_pinv @ x_full.T @ y
    adj = lmat @ beta
    vcov = lmat @ xtx_pinv @ lmat.T  # Var / mse
    d = np.diag(vcov)
    pair_var = pd.DataFrame(
        d[:, None] + d[None, :] - 2 * vcov, index=labels, columns=labels
    ).clip(lower=0.0)
    raw = df.groupby(treatment)[response].agg(["mean", "size"]).loc[labels]
    means = pd.DataFrame(
        {
            "treatment": labels,
            "adj_mean": adj,
            "se": np.sqrt(mse * d),
            "raw_mean": raw["mean"].to_numpy(),
            "n": raw["size"].to_numpy(),
        }
    ).sort_values("adj_mean", ascending=False, kind="stable").reset_index(drop=True)
    sed = np.sqrt(mse * pair_var)

    grand = float(y.mean())
    cv = float(np.sqrt(mse) / grand * 100) if grand != 0 else float("nan")
    notes = (
        "Model: intra-block analysis of an alpha lattice (fixed replicates and blocks\n"
        "  nested in replicates); treatments adjusted for blocks.\n"
        "  y_ijk = μ + ρ_i + β_j(i) + τ_k + ε_ijk\n"
        "Means are least-squares (adjusted) means; each pair has its own SED.\n"
        "Assumptions: independent plots, approximately normal residuals, homogeneous variances.\n"
        "Statistical significance does not imply agronomic importance.\n"
        "CV% = sqrt(MS_error) / grand_mean × 100."
    )
    res = ANOVAResult(
        table=table, design="Alpha lattice", n_obs=n, grand_mean=grand, cv_pct=cv,
        model_notes=notes, missing_reported=missing,
    )
    rel_eff = float(sse_rcbd / (n - r_rcbd) / mse * 100)
    return AlphaLatticeResult(
        anova=res, means=means, sed=sed, mse=float(mse), df_error=int(df_err),
        relative_efficiency=rel_eff, alpha=alpha, decimals=decimals, _pair_var=pair_var,
    )

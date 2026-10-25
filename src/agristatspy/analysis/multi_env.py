"""
Combined analysis over environments (locations, years, or location × year).

Model for a balanced RCBD in every environment
----------------------------------------------
    y_ijk = μ + E_i + R_j(i) + T_k + (TE)_ik + ε_ijk

Sources: Environment, Replication (Block) within environment, Treatment,
Treatment × Environment and pooled Error.  A CRD in each environment (no
``block``) drops the replication term.

Testing
-------
* ``environment_effect="random"`` (default; the usual choice for locations
  and years sampled from a wider population): treatment is tested against
  Treatment × Environment, the interaction against the pooled error.
* ``environment_effect="fixed"``: both are tested against the pooled error.

Environment is always tested against Replication within environment.  Before
pooling, error variances of the environments are compared with Bartlett's
test; if they differ markedly the combined analysis is unreliable and the
environments should be analysed separately (see ``per_env``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis._export import export_tables
from agristatspy.analysis.alpha_lattice import _add_letters
from agristatspy.analysis.anova import ANOVAResult, anova_crd, anova_rcbd
from agristatspy.analysis.comparison import ComparisonResult, compare_means
from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError

_SEP = " | "


@dataclass
class MultiEnvResult:
    """Result of :func:`multi_environment`.

    Attributes
    ----------
    per_env:
        ``{environment: ANOVAResult}`` — the separate analysis of each
        environment.
    homogeneity:
        Dict with Bartlett's ``chi2``, ``df``, ``p``, the largest-to-smallest
        error-MS ``ratio``, ``homogeneous`` (``p >= 0.05``) and an ``error_ms``
        table.
    combined:
        :class:`ANOVAResult` of the combined analysis.
    means:
        Treatment means over all environments: ``treatment``, ``mean``,
        ``sem`` (from the error term used for treatment), ``n``.
    mean_table:
        Treatment × environment table of means, with an ``Overall`` column.
    stats:
        ``mse``, ``df_error`` (pooled), ``mse_te``, ``df_te``, ``grand_mean``,
        ``cv_pct``, ``cd`` and ``alpha``.
    """

    design: str
    environment_effect: str
    per_env: dict[str, ANOVAResult]
    homogeneity: dict
    combined: ANOVAResult
    means: pd.DataFrame
    mean_table: pd.DataFrame
    stats: dict
    alpha: float = 0.05
    decimals: int = 2
    messages: list[str] = field(default_factory=list)
    comparisons: dict[str, ComparisonResult] = field(default_factory=dict)
    genetics: object | None = None
    _n_per_mean: int = 0

    def compare(self, method: str = "tukey", alpha: float | None = None) -> ComparisonResult:
        """Compare treatment means over environments.

        Uses the error term that tested treatment in the combined ANOVA:
        Treatment × Environment for random environments, pooled error for
        fixed.  ``.letters()`` gives compact letter groups.
        """
        s = self.stats
        if self.environment_effect == "random":
            mse, dfe, note = s["mse_te"], s["df_te"], "Treatment × Environment mean square"
        else:
            mse, dfe, note = s["mse"], s["df_error"], "pooled error over environments"
        return compare_means(
            self.means.set_index("treatment")["mean"], self._n_per_mean, mse, dfe,
            method=method, alpha=self.alpha if alpha is None else alpha, error_note=note,
        )

    def thesis_table(self, decimals: int | None = None) -> pd.DataFrame:
        d = self.decimals if decimals is None else decimals
        m = self.means
        out = pd.DataFrame(
            {
                "Treatment": m["treatment"].astype(str),
                "Mean ± SEM": [
                    f"{a:.{d}f} ± {b:.{d}f}" for a, b in zip(m["mean"], m["sem"], strict=True)
                ],
            }
        )
        return _add_letters(out, self.comparisons)

    def tables(self) -> dict[str, pd.DataFrame]:
        h = self.homogeneity
        out = {
            "Treatment means": self.thesis_table(),
            "Means by environment": self.mean_table.reset_index(),
            "Combined ANOVA": self.combined.table,
            "Error homogeneity": h["error_ms"],
        }
        for c in self.comparisons.values():
            out[f"Pairwise {c.method}"[:31]] = c.table
        if self.genetics is not None:
            out["Genetic parameters"] = self.genetics.table
        for env, a in self.per_env.items():
            out[f"ANOVA {env}"[:31]] = a.table
        return out

    def export(self, path: str) -> object:
        """Export to ``.xlsx``, ``.csv``, ``.md``, ``.html`` or ``.tex``."""
        return export_tables(self.tables(), path, "Combined analysis over environments", str(self))

    def __str__(self) -> str:
        h, s = self.homogeneity, self.stats
        lines = [
            f"Combined analysis over {len(self.per_env)} environments — {self.design} "
            f"({self.environment_effect} environments)",
            "",
            "1. Homogeneity of error variances (Bartlett)",
            f"   χ² = {h['chi2']:.3f}, df = {h['df']}, p = {h['p']:.4f}; "
            f"largest/smallest error MS = {h['ratio']:.2f}  →  "
            + ("homogeneous; pooling is reasonable" if h["homogeneous"]
               else "NOT homogeneous; interpret the combined analysis with caution"),
            "",
            h["error_ms"].to_string(index=False),
            "",
            "2. Combined ANOVA",
            f"   Grand mean = {s['grand_mean']:.4f}   CV = {s['cv_pct']:.2f}%   "
            f"CD({int(self.alpha * 100)}%) = {s['cd']:.4f}",
            "",
            self.combined.table.to_string(index=False, float_format=lambda x: f"{x:.4f}"),
            "",
            "3. Treatment means over environments"
            + (f" (letters: {', '.join(c.method for c in self.comparisons.values())})"
               if self.comparisons else ""),
            "",
            self.thesis_table().to_string(index=False),
        ]
        if self.messages:
            lines += ["", *[f"• {m}" for m in self.messages]]
        if self.genetics is not None:
            lines += ["", "4. " + str(self.genetics)]
        lines += ["", "Statistical significance does not imply agronomic importance."]
        return "\n".join(lines)


def _bartlett(mses: np.ndarray, dfs: np.ndarray) -> tuple[float, int, float]:
    k = len(mses)
    pooled = float((dfs * mses).sum() / dfs.sum())
    c = 1 + (np.sum(1 / dfs) - 1 / dfs.sum()) / (3 * (k - 1))
    chi2 = float((dfs.sum() * np.log(pooled) - np.sum(dfs * np.log(mses))) / c)
    return chi2, k - 1, float(stats.chi2.sf(chi2, k - 1))


def multi_environment(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    environment: str | list[str],
    block: str | None = None,
    environment_effect: str = "random",
    alpha: float = 0.05,
    decimals: int = 2,
) -> MultiEnvResult:
    """Combined ANOVA over locations, years, or location × year.

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per plot.
    response, treatment:
        Response and treatment columns.
    environment:
        Environment column (location or year), or a list of columns such as
        ``["location", "year"]`` that together define an environment.
    block:
        Block/replicate column (RCBD in each environment).  ``None`` means a
        CRD in each environment.  Labels may repeat across environments.
    environment_effect:
        ``"random"`` (default) tests treatment against Treatment ×
        Environment; ``"fixed"`` tests it against the pooled error.
    alpha:
        Significance level (used for the critical difference).
    decimals:
        Decimals in formatted tables.

    Returns
    -------
    MultiEnvResult

    Raises
    ------
    AnalysisError
        If columns are missing, fewer than 2 environments are present, or the
        layout is not complete and balanced (every treatment once per block in
        every environment, with the same number of replicates everywhere).

    Notes
    -----
    A significant Treatment × Environment interaction means treatment
    rankings change across environments; overall means then describe the
    average, not any one environment.  Genotype-by-environment stability
    methods (AMMI, GGE, regression) are not included.
    Statistical significance does not establish agronomic importance.
    """
    if environment_effect not in ("random", "fixed"):
        raise AnalysisError("environment_effect must be 'random' or 'fixed'.")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
        raise AnalysisError(f"alpha must be a number strictly between 0 and 1; got {alpha!r}.")
    env_cols = [environment] if isinstance(environment, str) else list(environment)
    cols = [response, treatment, *env_cols] + ([block] if block else [])
    _check_columns(data, cols, "multi_environment")
    _check_numeric(data, response, "multi_environment")

    df = data[cols].copy()
    n_before = len(df)
    df = df.dropna(subset=[response]).reset_index(drop=True)
    if n_before - len(df):
        raise AnalysisError(
            f"multi_environment: {n_before - len(df)} missing response(s); the combined "
            "analysis needs a complete, balanced layout."
        )
    df["__env__"] = df[env_cols].astype(str).agg(_SEP.join, axis=1)
    df["__trt__"] = df[treatment].astype(str)
    df["y"] = df[response].astype(float)
    envs = list(dict.fromkeys(df["__env__"]))
    trts = list(dict.fromkeys(df["__trt__"]))
    e, t = len(envs), len(trts)
    if e < 2:
        raise AnalysisError(f"multi_environment needs at least 2 environments; found {e}.")
    if t < 2:
        raise AnalysisError(f"multi_environment needs at least 2 treatments; found {t}.")

    if block:
        df["__blk__"] = df["__env__"] + "/" + df[block].astype(str)
        cell = df.groupby(["__env__", "__trt__", "__blk__"]).size()
        nblk = df.groupby("__env__")["__blk__"].nunique()
        if nblk.nunique() != 1:
            raise AnalysisError(
                f"multi_environment: environments have different numbers of blocks:\n{nblk}"
            )
        r = int(nblk.iloc[0])
        if len(cell) != e * t * r or cell.nunique() != 1:
            raise AnalysisError(
                "multi_environment: every treatment must appear exactly once in every "
                "block of every environment (complete, balanced layout)."
            )
    else:
        cell = df.groupby(["__env__", "__trt__"]).size()
        if len(cell) != e * t or cell.nunique() != 1:
            raise AnalysisError(
                "multi_environment: every environment × treatment needs the same number "
                "of replicates (complete, balanced layout)."
            )
        r = int(cell.iloc[0])
    if r < 2 and not block:
        raise AnalysisError("multi_environment needs at least 2 replicates per environment.")

    # ---- per-environment ANOVA and error homogeneity ----------------------
    per_env: dict[str, ANOVAResult] = {}
    for env in envs:
        sub = df[df["__env__"] == env]
        per_env[env] = (
            anova_rcbd(sub, "y", "__trt__", "__blk__") if block
            else anova_crd(sub, "y", "__trt__")
        )
    err = pd.DataFrame(
        {
            "environment": envs,
            "error_ms": [float(a.table.set_index("source").loc["Error", "ms"]) for a in per_env.values()],
            "error_df": [int(a.table.set_index("source").loc["Error", "df"]) for a in per_env.values()],
            "mean": [a.grand_mean for a in per_env.values()],
            "cv_pct": [a.cv_pct for a in per_env.values()],
        }
    )
    if (err["error_ms"] <= 0).any():
        raise AnalysisError("multi_environment: an environment has zero error variance.")
    chi2, hdf, hp = _bartlett(err["error_ms"].to_numpy(), err["error_df"].to_numpy(float))
    hom = {
        "chi2": chi2, "df": hdf, "p": hp, "homogeneous": hp >= 0.05,
        "ratio": float(err["error_ms"].max() / err["error_ms"].min()), "error_ms": err,
    }

    # ---- combined ANOVA (balanced sums of squares) -------------------------
    y = df["y"]
    gm = float(y.mean())
    env_m = df.groupby("__env__")["y"].mean()
    trt_m = df.groupby("__trt__")["y"].mean()
    cell_m = df.groupby(["__trt__", "__env__"])["y"].mean()
    ss_e = float(t * r * ((env_m - gm) ** 2).sum())
    ss_t = float(e * r * ((trt_m - gm) ** 2).sum())
    te = cell_m - trt_m.reindex(cell_m.index.get_level_values(0)).to_numpy() \
        - env_m.reindex(cell_m.index.get_level_values(1)).to_numpy() + gm
    ss_te = float(r * (te**2).sum())
    ss_tot = float(((y - gm) ** 2).sum())
    if block:
        blk_m = df.groupby(["__env__", "__blk__"])["y"].mean()
        ss_r = float(t * ((blk_m - env_m.reindex(blk_m.index.get_level_values(0)).to_numpy()) ** 2).sum())
        df_r = e * (r - 1)
    else:
        ss_r, df_r = 0.0, 0
    df_e, df_t, df_te = e - 1, t - 1, (t - 1) * (e - 1)
    df_err = e * (r - 1) * (t - 1) if block else e * t * (r - 1)
    ss_err = ss_tot - ss_e - ss_r - ss_t - ss_te
    mse = ss_err / df_err
    ms_te = ss_te / df_te
    nan = float("nan")

    def row(src, d, ss, den_ms=None, den_df=None):
        ms = ss / d
        f, p = (ms / den_ms, float(stats.f.sf(ms / den_ms, d, den_df))) if den_ms else (nan, nan)
        return {"source": src, "df": d, "ss": ss, "ms": ms, "f_stat": f, "p_value": p}

    rand = environment_effect == "random"
    rows = [row("Environment", df_e, ss_e, ss_r / df_r if block else None, df_r if block else None)]
    if block:
        rows.append(row("Block (Environment)", df_r, ss_r))
    rows += [
        row("Treatment", df_t, ss_t, ms_te if rand else mse, df_te if rand else df_err),
        row("Treatment × Environment", df_te, ss_te, mse, df_err),
        {"source": "Error", "df": df_err, "ss": ss_err, "ms": mse, "f_stat": nan, "p_value": nan},
        {"source": "Total", "df": len(df) - 1, "ss": ss_tot, "ms": nan, "f_stat": nan, "p_value": nan},
    ]
    cv = float(np.sqrt(mse) / gm * 100) if gm else nan
    design = "RCBD" if block else "CRD"
    notes = (
        f"Model: combined {design} analysis over environments "
        f"({'random' if rand else 'fixed'} environments).\n"
        + ("  Treatment is tested against Treatment × Environment.\n" if rand
           else "  Treatment is tested against the pooled error.\n")
        + "  Treatment × Environment is tested against the pooled error.\n"
        + ("  Environment is tested against Block (Environment).\n" if block else "")
        + "Assumptions: balanced layout, homogeneous error variances across environments,\n"
        "  independent plots, approximately normal residuals.\n"
        "Statistical significance does not imply agronomic importance.\n"
        "CV% = sqrt(MS_error) / grand_mean × 100."
    )
    combined = ANOVAResult(
        table=pd.DataFrame(rows), design=f"Combined {design}", n_obs=len(df),
        grand_mean=gm, cv_pct=cv, model_notes=notes,
    )

    n_per = r * e
    ms_for_t, df_for_t = (ms_te, df_te) if rand else (mse, df_err)
    sem = float(np.sqrt(ms_for_t / n_per))
    means = pd.DataFrame(
        {"treatment": trt_m.index, "mean": trt_m.to_numpy(), "sem": sem, "n": n_per}
    ).sort_values("mean", ascending=False, kind="stable").reset_index(drop=True)
    mt = df.pivot_table(index="__trt__", columns="__env__", values="y", aggfunc="mean")[envs]
    mt["Overall"] = trt_m
    mt.index.name = "treatment"
    mt = mt.loc[means["treatment"]]
    cd = float(stats.t.ppf(1 - alpha / 2, df_for_t) * np.sqrt(2 * ms_for_t / n_per))

    msgs = []
    if not hom["homogeneous"]:
        msgs.append(
            f"Error variances differ across environments (Bartlett p = {hp:.4f}); "
            "consider analysing environments separately (see .per_env)."
        )
    p_te = rows[-3]["p_value"]
    if p_te < alpha:
        msgs.append(
            "Treatment × Environment is significant: treatment performance depends on the "
            "environment, so overall means should be interpreted with care."
        )
    return MultiEnvResult(
        design=design, environment_effect=environment_effect, per_env=per_env,
        homogeneity=hom, combined=combined, means=means, mean_table=mt,
        stats={"mse": float(mse), "df_error": int(df_err), "mse_te": float(ms_te),
               "df_te": int(df_te), "grand_mean": gm, "cv_pct": cv, "cd": cd, "alpha": alpha},
        alpha=alpha, decimals=decimals, messages=msgs, _n_per_mean=n_per,
    )

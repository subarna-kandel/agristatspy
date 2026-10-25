"""
Genetic variability parameters from replicated trials.

Variance components (expected mean squares)
-------------------------------------------
Single environment, r replications:
    σ²g = (MSg − MSe) / r        σ²e = MSe        σ²p = σ²g + σ²e
Over e environments (random), r replications:
    σ²ge = (MSge − MSe) / r      σ²g = (MSg − MSge) / (r·e)
    σ²p = σ²g + σ²ge / e + σ²e / (r·e)

Broad-sense heritability h² = σ²g / σ²p (on a genotype-mean basis);
GCV = 100·σg / mean; PCV = 100·σp / mean; ECV = 100·σe / mean;
genetic advance GA = k·σp·h² (k = 2.06 for 5 % selection intensity) and
GAM = 100·GA / mean.  A negative variance component is set to zero.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from agristatspy.analysis._export import export_tables
from agristatspy.analysis.anova import anova_crd, anova_rcbd
from agristatspy.analysis.multi_env import multi_environment
from agristatspy.analysis.summary import _check_columns, _check_numeric
from agristatspy.exceptions import AnalysisError


def _band(v: float, low: float, high: float) -> str:
    if np.isnan(v):
        return ""
    return "low" if v < low else "moderate" if v < high else "high"


@dataclass
class GeneticParameters:
    """Result of :func:`genetic_parameters`.

    ``table`` has one row per trait with the mean squares, variance
    components, GCV, PCV, ECV, ``h2_bs`` (broad-sense heritability, %),
    ``ga``, ``gam`` and the usual qualitative classes (h²: Robinson 1966;
    GAM: Johnson et al. 1955; GCV/PCV: Deshmukh et al. 1986).
    """

    table: pd.DataFrame
    k: float
    notes: list[str]

    def tables(self) -> dict[str, pd.DataFrame]:
        return {"Genetic parameters": self.table}

    def export(self, path: str) -> object:
        return export_tables(self.tables(), path, "Genetic parameters", str(self))

    def __str__(self) -> str:
        show = self.table[
            ["trait", "mean", "gcv", "pcv", "h2_bs", "ga", "gam", "h2_class", "gam_class"]
        ]
        lines = ["Genetic parameters", "", show.to_string(index=False, float_format=lambda x: f"{x:.3f}")]
        lines += ["", f"GA uses selection differential k = {self.k}.", *self.notes,
                  "Heritability is specific to this population and environment(s)."]
        return "\n".join(lines)


def genetic_parameters(
    data: pd.DataFrame,
    traits: list[str],
    genotype: str,
    replication: str | None = None,
    environment: str | list[str] | None = None,
    k: float = 2.06,
) -> GeneticParameters:
    """Genotypic and phenotypic variance, GCV, PCV, heritability, genetic advance.

    Parameters
    ----------
    data:
        Tidy DataFrame, one row per plot.
    traits:
        Trait columns.
    genotype:
        Genotype (treatment) column; genotypes are treated as random.
    replication:
        Replication/block column (RCBD).  ``None`` analyses a CRD.
    environment:
        Optional environment column(s).  When given, the combined analysis
        over environments is used and the genotype × environment variance is
        separated.
    k:
        Standardized selection differential (2.06 = top 5 %, 1.76 = 10 %,
        1.40 = 20 %).

    Returns
    -------
    GeneticParameters

    Raises
    ------
    AnalysisError
        On missing columns, fewer than 2 genotypes, or incomplete layouts.

    Notes
    -----
    Components are method-of-moments estimates from a balanced fixed-effects
    ANOVA table; they are not REML estimates.  Heritability estimates depend
    on the population and environments sampled.
    """
    if isinstance(traits, str) or not len(traits):
        raise AnalysisError("genetic_parameters: 'traits' must be a non-empty list of columns.")
    cols = [genotype, *traits] + ([replication] if replication else [])
    _check_columns(data, cols, "genetic_parameters")
    for tr in traits:
        _check_numeric(data, tr, "genetic_parameters")
    if data[genotype].nunique() < 2:
        raise AnalysisError("genetic_parameters needs at least 2 genotypes.")

    rows, notes = [], []
    for tr in traits:
        if environment is not None:
            res = multi_environment(
                data, tr, genotype, environment, block=replication,
                environment_effect="random",
            )
            t = res.combined.table.set_index("source")
            e = len(res.per_env)
            r = res._n_per_mean // e
            ms_g, ms_ge, ms_e = (
                float(t.loc["Treatment", "ms"]),
                float(t.loc["Treatment × Environment", "ms"]),
                float(t.loc["Error", "ms"]),
            )
            raw_ge = (ms_ge - ms_e) / r
            raw_g = (ms_g - ms_ge) / (r * e)
            s2ge, s2g, s2e = max(raw_ge, 0.0), max(raw_g, 0.0), ms_e
            s2p = s2g + s2ge / e + s2e / (r * e)
            mean = res.combined.grand_mean
            extra = {"ms_ge": ms_ge, "s2ge": s2ge}
            neg = raw_g < 0 or raw_ge < 0
        else:
            sub = data[[genotype, tr] + ([replication] if replication else [])]
            a = anova_rcbd(sub, tr, genotype, replication) if replication else anova_crd(sub, tr, genotype)
            t = a.table.set_index("source")
            r = int(sub.dropna(subset=[tr]).groupby(genotype).size().iloc[0])
            ms_g, ms_e = float(t.loc[genotype if genotype in t.index else "Treatment", "ms"]), float(t.loc["Error", "ms"])
            raw_g = (ms_g - ms_e) / r
            s2g, s2e = max(raw_g, 0.0), ms_e
            s2p = s2g + s2e
            mean = a.grand_mean
            extra = {}
            neg = raw_g < 0
        if neg:
            notes.append(f"{tr}: a negative variance component was set to zero.")
        sg, sp, se = np.sqrt(s2g), np.sqrt(s2p), np.sqrt(s2e)
        h2 = s2g / s2p * 100 if s2p > 0 else np.nan
        ga = k * sp * h2 / 100 if s2p > 0 else np.nan
        scale = 100 / mean if mean else np.nan
        gam = ga * scale
        rows.append(
            {
                "trait": tr, "mean": mean, "ms_genotype": ms_g, "ms_error": ms_e, **extra,
                "var_g": s2g, "var_e": s2e, "var_p": s2p,
                "gcv": sg * scale, "pcv": sp * scale, "ecv": se * scale,
                "h2_bs": h2, "ga": ga, "gam": gam,
                "h2_class": _band(h2, 30, 60), "gam_class": _band(gam, 10, 20),
                "gcv_class": _band(sg * scale, 10, 20), "pcv_class": _band(sp * scale, 10, 20),
            }
        )
    return GeneticParameters(table=pd.DataFrame(rows), k=k, notes=notes)

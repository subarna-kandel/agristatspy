"""
Correlation and principal component analysis across traits.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from agristatspy.analysis._export import export_tables
from agristatspy.analysis.summary import _check_columns
from agristatspy.exceptions import AnalysisError


def _trait_frame(
    data: pd.DataFrame, traits: list[str] | None, by: str | None, context: str
) -> pd.DataFrame:
    if traits is None:
        traits = [
            c for c in data.columns
            if pd.api.types.is_numeric_dtype(data[c]) and c != by
        ]
    if isinstance(traits, str) or len(traits) < 2:
        raise AnalysisError(f"{context}: need at least 2 trait columns; got {traits!r}.")
    _check_columns(data, list(traits) + ([by] if by else []), context)
    for c in traits:
        if not pd.api.types.is_numeric_dtype(data[c]):
            raise AnalysisError(f"{context}: trait '{c}' must be numeric.")
    df = data[list(traits) + ([by] if by else [])]
    if by:
        df = df.groupby(by)[list(traits)].mean()
    return df[list(traits)]


@dataclass
class CorrelationResult:
    """Pairwise correlations between traits.

    ``r`` and ``p`` are trait × trait DataFrames, ``n`` the number of
    complete pairs.  P-values are unadjusted.
    """

    r: pd.DataFrame
    p: pd.DataFrame
    n: pd.DataFrame
    method: str

    def stars(self) -> pd.DataFrame:
        """``r`` formatted with significance stars (* 0.05, ** 0.01, *** 0.001)."""
        def f(r: float, p: float) -> str:
            s = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            return f"{r:.3f}{s}"
        out = self.r.copy().astype(object)
        for i in self.r.index:
            for j in self.r.columns:
                out.loc[i, j] = "1" if i == j else f(self.r.loc[i, j], self.p.loc[i, j])
        return out

    def table(self) -> pd.DataFrame:
        """Long format: ``trait1``, ``trait2``, ``r``, ``p``, ``n``."""
        rows = [
            {"trait1": a, "trait2": b, "r": self.r.loc[a, b], "p": self.p.loc[a, b], "n": self.n.loc[a, b]}
            for i, a in enumerate(self.r.index) for b in self.r.columns[i + 1:]
        ]
        return pd.DataFrame(rows)

    def plot(self, **kwargs: object):
        """Heat map of ``r`` with significance stars; returns ``(fig, ax)``."""
        from agristatspy.visualization.plots import plot_correlation

        return plot_correlation(self, **kwargs)

    def export(self, path: str) -> object:
        tabs = {"Correlation": self.stars().reset_index(names="trait"),
                "Pairs": self.table()}
        return export_tables(tabs, path, f"{self.method.title()} correlation", str(self))

    def __str__(self) -> str:
        return (
            f"{self.method.title()} correlation (* p<0.05, ** p<0.01, *** p<0.001)\n\n"
            + self.stars().to_string()
        )


def correlation(
    data: pd.DataFrame,
    traits: list[str] | None = None,
    method: str = "pearson",
    by: str | None = None,
) -> CorrelationResult:
    """Correlation matrix of traits with p-values.

    Parameters
    ----------
    data:
        DataFrame with one column per trait.
    traits:
        Trait columns (default: all numeric columns).
    method:
        ``"pearson"`` (default) or ``"spearman"``.
    by:
        Optional column (e.g. treatment); traits are averaged within each
        level first, so the correlation is between treatment means.

    Raises
    ------
    AnalysisError
        On unknown method, non-numeric traits, or fewer than 2 traits.

    Notes
    -----
    Correlation is not causation, and a significant coefficient between a few
    treatment means can be driven by one point.  Pairwise-complete data are
    used for each pair.  P-values are not adjusted for multiple testing.
    """
    if method not in ("pearson", "spearman"):
        raise AnalysisError("method must be 'pearson' or 'spearman'.")
    df = _trait_frame(data, traits, by, "correlation")
    fn = stats.pearsonr if method == "pearson" else stats.spearmanr
    names = list(df.columns)
    r = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    p = pd.DataFrame(np.zeros((len(names),) * 2), index=names, columns=names)
    n = pd.DataFrame(0, index=names, columns=names)
    for i, a in enumerate(names):
        n.loc[a, a] = int(df[a].notna().sum())
        for b in names[i + 1:]:
            pair = df[[a, b]].dropna()
            n.loc[a, b] = n.loc[b, a] = len(pair)
            if len(pair) < 3 or pair[a].nunique() < 2 or pair[b].nunique() < 2:
                rv, pv = np.nan, np.nan
            else:
                rv, pv = (float(x) for x in fn(pair[a], pair[b]))
            r.loc[a, b] = r.loc[b, a] = rv
            p.loc[a, b] = p.loc[b, a] = pv
    return CorrelationResult(r=r, p=p, n=n, method=method)


@dataclass
class PCAResult:
    """Principal component analysis.

    Attributes
    ----------
    eigenvalues:
        ``component``, ``eigenvalue``, ``variance_pct``, ``cumulative_pct``.
    loadings:
        Trait × component correlations between traits and components
        (eigenvector × √eigenvalue) when ``scale=True``; otherwise the
        eigenvectors scaled the same way.
    eigenvectors:
        Trait × component unit eigenvectors.
    scores:
        Observation (or group-mean) × component scores.
    scale:
        Whether traits were standardized (PCA of the correlation matrix).
    """

    eigenvalues: pd.DataFrame
    loadings: pd.DataFrame
    eigenvectors: pd.DataFrame
    scores: pd.DataFrame
    scale: bool
    n_dropped: int = 0

    def plot(self, **kwargs: object):
        """Scree plot and biplot of the first two components; ``(fig, axes)``."""
        from agristatspy.visualization.plots import plot_pca

        return plot_pca(self, **kwargs)

    def export(self, path: str) -> object:
        tabs = {
            "Eigenvalues": self.eigenvalues,
            "Loadings": self.loadings.reset_index(names="trait"),
            "Scores": self.scores.reset_index(),
        }
        return export_tables(tabs, path, "Principal component analysis", str(self))

    def __str__(self) -> str:
        k = min(4, self.loadings.shape[1])
        lines = [
            "Principal component analysis"
            + (" (standardized traits; correlation matrix)" if self.scale else " (covariance matrix)"),
            "",
            self.eigenvalues.to_string(index=False, float_format=lambda x: f"{x:.3f}"),
            "",
            f"Loadings (first {k} components):",
            self.loadings.iloc[:, :k].to_string(float_format=lambda x: f"{x:.3f}"),
        ]
        if self.n_dropped:
            lines += ["", f"{self.n_dropped} row(s) with missing values were dropped."]
        return "\n".join(lines)


def pca(
    data: pd.DataFrame,
    traits: list[str] | None = None,
    by: str | None = None,
    scale: bool = True,
) -> PCAResult:
    """Principal component analysis of traits.

    Parameters
    ----------
    data:
        DataFrame with one column per trait.
    traits:
        Trait columns (default: all numeric columns).
    by:
        Optional column (e.g. genotype); traits are averaged within each level
        and the levels become the observations.
    scale:
        Standardize traits to unit variance (default), i.e. PCA of the
        correlation matrix.  Use ``False`` only when traits share a unit.

    Raises
    ------
    AnalysisError
        If fewer than 3 complete observations or a trait has no variation.

    Notes
    -----
    Rows with missing values are dropped.  The first component's sign is
    arbitrary; it is fixed so the largest loading in each component is
    positive.
    """
    df = _trait_frame(data, traits, by, "pca")
    n0 = len(df)
    df = df.dropna()
    if len(df) < 3:
        raise AnalysisError("pca needs at least 3 complete observations.")
    if (df.std(ddof=1) == 0).any():
        raise AnalysisError("pca: a trait has no variation.")
    x = df - df.mean()
    if scale:
        x = x / df.std(ddof=1)
    cov = np.cov(x.to_numpy(), rowvar=False, ddof=1)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = np.clip(vals[order], 0, None), vecs[:, order]
    for j in range(vecs.shape[1]):
        if vecs[np.abs(vecs[:, j]).argmax(), j] < 0:
            vecs[:, j] = -vecs[:, j]
    ncomp = min(len(df) - 1, df.shape[1])
    vals, vecs = vals[:ncomp], vecs[:, :ncomp]
    names = [f"PC{i + 1}" for i in range(ncomp)]
    total = float(np.trace(cov))
    eig = pd.DataFrame(
        {
            "component": names, "eigenvalue": vals,
            "variance_pct": vals / total * 100,
            "cumulative_pct": np.cumsum(vals) / total * 100,
        }
    )
    return PCAResult(
        eigenvalues=eig,
        loadings=pd.DataFrame(vecs * np.sqrt(vals), index=df.columns, columns=names),
        eigenvectors=pd.DataFrame(vecs, index=df.columns, columns=names),
        scores=pd.DataFrame(x.to_numpy() @ vecs, index=df.index, columns=names),
        scale=scale, n_dropped=n0 - len(df),
    )

"""
Matplotlib plots: field layouts, treatment means, boxplots, residual diagnostics.

All functions return ``(fig, ax)`` (or ``(fig, axes)``) and never call
``plt.show()``; save with ``fig.savefig("name.png", dpi=300)`` or ``.svg``.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from agristatspy.exceptions import AnalysisError


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _colors(labels: list[str]) -> dict[str, tuple]:
    plt = _plt()
    cmap = plt.get_cmap("tab20")
    return {lab: cmap(i % 20) for i, lab in enumerate(labels)}


def plot_layout(design: object, figsize: tuple[float, float] | None = None):
    """Draw a field layout: one cell per plot, coloured by treatment.

    RCBD, factorial-RCBD and split-plot designs are drawn with one band per
    block (with a gap between blocks); Latin squares as a row × column grid;
    CRD as a grid using the design's ``row``/``column`` (or a grid as wide as
    the number of treatments).  The picture is a view of the table and never
    changes the allocation.

    Returns
    -------
    (Figure, Axes)
    """
    plt = _plt()
    meta = design.metadata  # type: ignore[attr-defined]
    tbl = design.table  # type: ignore[attr-defined]
    base = meta.extra.get("base_design", meta.design_type)

    positions: list[tuple[float, float]] = []  # (x, y) of each cell's lower-left
    block_boxes: list[tuple[float, float, float, float, str]] = []
    gap = 0.6

    if base == "LSQ":
        for _, r in tbl.iterrows():
            positions.append((r["column"] - 1, -(r["row"] - 1) - 1))
        width, height = meta.n_treatments, meta.n_treatments
    elif base == "CRD":
        if "row" in tbl.columns:
            for _, r in tbl.iterrows():
                positions.append((r["column"] - 1, -(r["row"] - 1) - 1))
            width, height = int(tbl["column"].max()), int(tbl["row"].max())
        else:
            ncol = meta.n_treatments
            for i in range(len(tbl)):
                positions.append((i % ncol, -(i // ncol) - 1))
            width, height = ncol, math.ceil(len(tbl) / ncol)
    elif base == "AlphaLattice":
        k = meta.extra["block_size"]
        y_cursor = 0.0
        width = k
        for rep in sorted(tbl["replicate"].unique()):
            rdf = tbl[tbl["replicate"] == rep]
            top = y_cursor
            for blk in sorted(rdf["block"].unique()):
                bdf = rdf[rdf["block"] == blk].sort_values("position_in_block")
                for i in range(len(bdf)):
                    positions.append((i, top - 1))
                block_boxes.append((-0.05, top - 1.05, k + 0.1, 1.1, f"Rep {rep} · Block {blk}"))
                top -= 1.5
            y_cursor = top - gap
        height = -y_cursor
        tbl = tbl.sort_values(["replicate", "block", "position_in_block"])
    else:  # RCBD / SplitPlot: one band per block
        if base == "SplitPlot":
            per_block = len(meta.extra["main_levels"]) * len(meta.extra["sub_levels"])
            ncol = len(meta.extra["sub_levels"])
        else:
            per_block = meta.n_treatments
            ncol = per_block
        y_cursor = 0.0
        width = ncol
        for b in sorted(tbl["block"].unique()):
            bdf = tbl[tbl["block"] == b]
            nrows = math.ceil(per_block / ncol)
            top = y_cursor
            for i in range(len(bdf)):
                positions.append((i % ncol, top - (i // ncol) - 1))
            block_boxes.append((-0.05, top - nrows - 0.05, ncol + 0.1, nrows + 0.1, f"Block {b}"))
            y_cursor = top - nrows - gap
        height = -y_cursor

    labels = sorted(tbl["treatment"].unique())
    colors = _colors(labels)
    if figsize is None:
        figsize = (max(5.0, width * 1.5), max(3.0, height * 1.1 + 1))
    fig, ax = plt.subplots(figsize=figsize)
    for (x, y), (_, r) in zip(positions, tbl.iterrows(), strict=True):
        ax.add_patch(
            plt.Rectangle((x, y), 1, 1, facecolor=colors[r["treatment"]], edgecolor="black", lw=1)
        )
        ax.text(
            x + 0.5, y + 0.5, f"{r['plot_id']}\n{r['treatment']}",
            ha="center", va="center", fontsize=7,
        )
    for x, y, w, h, name in block_boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=False, edgecolor="black", lw=2.2))
        ax.text(x, y + h + 0.05, name, fontsize=8, fontweight="bold", va="bottom")

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    ax.set_xlim(min(xs) - 0.3, max(xs) + 1.3)
    ax.set_ylim(min(ys) - 0.3, max(ys) + 1.7 if block_boxes else max(ys) + 1.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        f"{meta.design_type}: {meta.n_treatments} treatments, {meta.total_plots} plots"
    )
    return fig, ax


def _group_stats(data: pd.DataFrame, response: str, treatment: str, error: str):
    from scipy import stats

    for c in (response, treatment):
        if c not in data.columns:
            raise AnalysisError(f"Column {c!r} not found in data.")
    g = data.dropna(subset=[response]).groupby(treatment)[response]
    mean, sd, n = g.mean(), g.std(ddof=1), g.size()
    se = sd / n.pow(0.5)
    if error == "se":
        err = se
    elif error == "sd":
        err = sd
    elif error == "ci95":
        err = se * pd.Series({k: stats.t.ppf(0.975, max(v - 1, 1)) for k, v in n.items()})
    else:
        raise AnalysisError("error must be 'se', 'sd', or 'ci95'.")
    return mean, err.fillna(0.0)


def plot_means(
    data: pd.DataFrame,
    response: str,
    treatment: str,
    error: str = "se",
    letters: pd.DataFrame | None = None,
):
    """Plot treatment means with error bars.

    Parameters
    ----------
    error:
        ``"se"`` (standard error, default), ``"sd"``, or ``"ci95"`` (per-group
        t-based 95% interval).  The choice is stated in the y-axis label.
    letters:
        Optional output of ``ComparisonResult.letters()`` to annotate bars.

    Returns
    -------
    (Figure, Axes)
    """
    plt = _plt()
    mean, err = _group_stats(data, response, treatment, error)
    fig, ax = plt.subplots(figsize=(max(4.0, 1.1 * len(mean) + 2), 4))
    colors = _colors(list(mean.index))
    ax.bar(
        range(len(mean)), mean.values, yerr=err.values, capsize=4,
        color=[colors[k] for k in mean.index], edgecolor="black",
    )
    ax.set_xticks(range(len(mean)))
    ax.set_xticklabels([str(k) for k in mean.index], rotation=30, ha="right")
    names = {"se": "± SE", "sd": "± SD", "ci95": "± 95% CI"}
    ax.set_ylabel(f"{response} (mean {names[error]})")
    ax.set_xlabel(treatment)
    if letters is not None:
        lmap = dict(zip(letters["treatment"], letters["group"], strict=True))
        for i, k in enumerate(mean.index):
            ax.text(i, mean[k] + err[k], f" {lmap.get(k, '')}", ha="center", va="bottom")
    ax.set_title(f"Treatment means of {response}")
    fig.tight_layout()
    return fig, ax


def plot_boxplot(data: pd.DataFrame, response: str, treatment: str):
    """Boxplot of the response by treatment, with individual points overlaid.

    Returns
    -------
    (Figure, Axes)
    """
    plt = _plt()
    for c in (response, treatment):
        if c not in data.columns:
            raise AnalysisError(f"Column {c!r} not found in data.")
    df = data.dropna(subset=[response])
    keys = sorted(df[treatment].unique(), key=str)
    vals = [df.loc[df[treatment] == k, response].astype(float).values for k in keys]
    fig, ax = plt.subplots(figsize=(max(4.0, 1.1 * len(keys) + 2), 4))
    ax.boxplot(vals, tick_labels=[str(k) for k in keys], showfliers=True)
    for i, v in enumerate(vals, start=1):
        ax.scatter([i] * len(v), v, s=14, color="black", alpha=0.6, zorder=3)
    ax.set_ylabel(response)
    ax.set_xlabel(treatment)
    ax.set_title(f"{response} by {treatment}")
    fig.tight_layout()
    return fig, ax


def plot_interaction(data: pd.DataFrame, response: str, x: str, trace: str):
    """Interaction plot: mean response across *x*, one line per *trace* level.

    Roughly parallel lines suggest little interaction; crossing or diverging
    lines suggest interaction.  This is a visual aid, not a test.

    Returns
    -------
    (Figure, Axes)
    """
    plt = _plt()
    for c in (response, x, trace):
        if c not in data.columns:
            raise AnalysisError(f"Column {c!r} not found in data.")
    m = data.dropna(subset=[response]).groupby([x, trace])[response].mean().unstack(trace)
    fig, ax = plt.subplots(figsize=(6, 4))
    for col in m.columns:
        ax.plot([str(i) for i in m.index], m[col].values, marker="o", label=str(col))
    ax.set_xlabel(x)
    ax.set_ylabel(f"mean {response}")
    ax.legend(title=trace)
    ax.set_title(f"Interaction plot: {x} × {trace}")
    fig.tight_layout()
    return fig, ax


def plot_residuals(diagnostics: object):
    """Residuals-vs-fitted and normal Q–Q plot for a ``DiagnosticsResult``.

    Returns
    -------
    (Figure, ndarray of Axes)
    """
    from scipy import stats

    plt = _plt()
    res = diagnostics.residuals  # type: ignore[attr-defined]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(res["fitted"], res["residual"], color="black", s=18)
    axes[0].axhline(0, color="grey", lw=1)
    axes[0].set_xlabel("Fitted value")
    axes[0].set_ylabel("Residual")
    axes[0].set_title("Residuals vs fitted")
    stats.probplot(res["residual"], dist="norm", plot=axes[1])
    axes[1].set_title("Normal Q–Q plot")
    fig.tight_layout()
    return fig, axes


def plot_correlation(result: object, figsize: tuple[float, float] | None = None):
    """Heat map of a :class:`~agristatspy.analysis.multivariate.CorrelationResult`.

    Cells show ``r`` with significance stars (* 0.05, ** 0.01, *** 0.001).
    Returns ``(Figure, Axes)``.
    """
    plt = _plt()
    r = result.r  # type: ignore[attr-defined]
    k = len(r)
    fig, ax = plt.subplots(figsize=figsize or (max(4.0, 0.9 * k + 2), max(3.5, 0.9 * k + 1)))
    im = ax.imshow(r.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(k), r.columns, rotation=45, ha="right")
    ax.set_yticks(range(k), r.index)
    stars = result.stars()  # type: ignore[attr-defined]
    for i in range(k):
        for j in range(k):
            ax.text(j, i, stars.iloc[i, j], ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="r")
    ax.set_title(f"{result.method.title()} correlation")  # type: ignore[attr-defined]
    fig.tight_layout()
    return fig, ax


def plot_pca(result: object, figsize: tuple[float, float] | None = None):
    """Scree plot and PC1-PC2 biplot of a :class:`PCAResult`.

    Points are observations (or group means); arrows are trait loadings.
    Returns ``(Figure, (ax_scree, ax_biplot))``.
    """
    plt = _plt()
    eig, load, sc = result.eigenvalues, result.loadings, result.scores  # type: ignore[attr-defined]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=figsize or (11, 4.5))
    a1.bar(eig["component"], eig["variance_pct"], color="#4c78a8")
    a1.plot(eig["component"], eig["cumulative_pct"], "o-", color="black")
    a1.set_ylabel("Variance explained (%)")
    a1.set_title("Scree plot")
    if sc.shape[1] >= 2:
        a2.scatter(sc["PC1"], sc["PC2"], color="#4c78a8")
        for lab, (x, y) in sc[["PC1", "PC2"]].iterrows():
            a2.annotate(str(lab), (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
        scale = float(np.abs(sc[["PC1", "PC2"]].to_numpy()).max() / max(np.abs(load[["PC1", "PC2"]].to_numpy()).max(), 1e-12))
        for lab, (x, y) in load[["PC1", "PC2"]].iterrows():
            a2.arrow(0, 0, x * scale * 0.9, y * scale * 0.9, color="#e45756", head_width=0.05 * scale * 0.3)
            a2.text(x * scale, y * scale, str(lab), color="#e45756", fontsize=8)
        a2.axhline(0, color="grey", lw=0.5)
        a2.axvline(0, color="grey", lw=0.5)
        a2.set_xlabel(f"PC1 ({eig['variance_pct'].iloc[0]:.1f}%)")
        a2.set_ylabel(f"PC2 ({eig['variance_pct'].iloc[1]:.1f}%)")
    a2.set_title("Biplot")
    fig.tight_layout()
    return fig, (a1, a2)

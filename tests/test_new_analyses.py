import numpy as np
import pandas as pd
import pytest
from scipy import stats

import agristatspy as ag
from agristatspy.exceptions import AnalysisError, DesignError


def _lattice_data(seed=0, t=12, k=4, r=3):
    d = ag.design.alpha_lattice(n_treatments=t, block_size=k, replications=r, seed=1).table
    rng = np.random.default_rng(seed)
    eff = {x: rng.normal(0, 3) for x in d["treatment"].unique()}
    beff = {(a, b): rng.normal(0, 2) for a in range(1, r + 1) for b in range(1, t // k + 1)}
    d["y"] = [50 + eff[x] + beff[(a, b)] + rng.normal(0, 1)
              for x, a, b in zip(d["treatment"], d["replicate"], d["block"], strict=True)]
    return d


def _env_data(seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for loc in "AB":
        for yr in (2023, 2024):
            for b in (1, 2, 3):
                be = rng.normal(0, 1)
                for i, t in enumerate(["T1", "T2", "T3", "T4"]):
                    rows.append({"loc": loc, "year": yr, "block": b, "trt": t,
                                 "y": 20 + 2 * i + be + rng.normal(0, 1), "z": rng.normal(10, 2)})
    return pd.DataFrame(rows)


# ---- alpha lattice -------------------------------------------------------

def test_alpha_lattice_design_valid_and_seeded():
    a = ag.design.alpha_lattice(n_treatments=12, block_size=4, replications=3, seed=5)
    b = ag.design.alpha_lattice(n_treatments=12, block_size=4, replications=3, seed=5)
    assert a.validate()
    assert a.table.equals(b.table)
    assert (a.table.groupby(["replicate", "treatment"]).size() == 1).all()
    assert "REPLICATE 1" in a.text_layout()
    fig, _ = a.plot_layout()
    assert fig is not None


@pytest.mark.parametrize("kw", [
    {"n_treatments": 10, "block_size": 4}, {"n_treatments": 12, "block_size": 1},
    {"n_treatments": 12, "block_size": 4, "replications": 4},
])
def test_alpha_lattice_design_errors(kw):
    with pytest.raises(DesignError):
        ag.design.alpha_lattice(**kw)


def test_alpha_lattice_anova_matches_ols():
    import statsmodels.formula.api as smf

    d = _lattice_data()
    res = ag.anova_alpha_lattice(d, "y", "treatment", "replicate", "block")
    d["blk"] = d["replicate"].astype(str) + "_" + d["block"].astype(str)
    m = smf.ols("y ~ C(blk) + C(treatment)", d).fit()
    assert res.mse == pytest.approx(m.mse_resid)
    assert res.df_error == m.df_resid
    t = res.anova.table.set_index("source")
    assert t.loc["Treatment (adj.)", "df"] == 11
    assert t["df"].iloc[:-1].sum() == t.loc["Total", "df"]
    assert res.relative_efficiency > 100
    c = res.compare("tukey")
    assert len(c.table) == 66
    assert set(c.letters().columns) == {"treatment", "mean", "group"}
    assert res.sed.shape == (12, 12)


def test_alpha_lattice_missing_value_ok():
    d = _lattice_data()
    d.loc[3, "y"] = np.nan
    res = ag.anova_alpha_lattice(d, "y", "treatment", "replicate", "block")
    assert res.anova.missing_reported == 1
    assert res.df_error == 15


# ---- multi environment ---------------------------------------------------

def test_multi_env_ss_and_tests():
    df = _env_data()
    res = ag.multi_environment(df, "y", "trt", ["loc", "year"], block="block")
    t = res.combined.table.set_index("source")
    assert t.loc["Treatment", "df"] == 3
    assert t.loc["Treatment × Environment", "df"] == 9
    assert t.loc["Error", "df"] == 4 * 2 * 3
    assert t["ss"].iloc[:-1].sum() == pytest.approx(t.loc["Total", "ss"])
    assert t.loc["Treatment", "f_stat"] == pytest.approx(
        t.loc["Treatment", "ms"] / t.loc["Treatment × Environment", "ms"])
    fixed = ag.multi_environment(df, "y", "trt", ["loc", "year"], block="block",
                                 environment_effect="fixed").combined.table.set_index("source")
    assert fixed.loc["Treatment", "f_stat"] == pytest.approx(
        fixed.loc["Treatment", "ms"] / fixed.loc["Error", "ms"])
    assert 0 <= res.homogeneity["p"] <= 1
    assert len(res.per_env) == 4
    assert list(res.compare("lsd").letters().columns) == ["treatment", "mean", "group"]


def test_multi_env_single_column_crd_and_errors():
    df = _env_data()
    one = ag.multi_environment(df[df["loc"] == "A"], "y", "trt", "year", block="block")
    assert len(one.per_env) == 2
    with pytest.raises(AnalysisError):
        ag.multi_environment(df[df["year"] == 2023].query("loc == 'A'"), "y", "trt", "year", block="block")
    with pytest.raises(AnalysisError):
        ag.multi_environment(df.iloc[:-1], "y", "trt", ["loc", "year"], block="block")


def test_multi_env_export(tmp_path):
    res = ag.multi_environment(_env_data(), "y", "trt", ["loc", "year"], block="block")
    for ext in ("xlsx", "csv", "md", "html", "tex"):
        assert res.export(tmp_path / f"r.{ext}").exists()


# ---- multi-trait / genetics ----------------------------------------------

def test_analyze_traits():
    df = _env_data()
    df = df[(df["loc"] == "A") & (df["year"] == 2023)]
    rep = ag.analyze_traits(df, ["y", "z"], "rcbd", treatment="trt", block="block")
    assert list(rep.summary["trait"]) == ["y", "z"]
    assert "p_Treatment" in rep.summary
    auto = ag.analyze_traits(df.drop(columns=["year"]), None, "rcbd", treatment="trt", block="block")
    assert set(auto.reports) == {"y", "z"}
    with pytest.raises(AnalysisError):
        ag.analyze_traits(df, ["nope"], "rcbd", treatment="trt", block="block")


def test_genetic_parameters_formulas():
    df = _env_data()
    df = df[(df["loc"] == "A") & (df["year"] == 2023)]
    gp = ag.genetic_parameters(df, ["y"], "trt", replication="block").table.iloc[0]
    a = ag.anova_rcbd(df, "y", "trt", "block").table.set_index("source")
    msg, mse = a.loc["Treatment", "ms"], a.loc["Error", "ms"]
    s2g = (msg - mse) / 3
    assert gp["var_g"] == pytest.approx(s2g)
    assert gp["h2_bs"] == pytest.approx(100 * s2g / (s2g + mse))
    assert gp["gcv"] == pytest.approx(100 * np.sqrt(s2g) / df["y"].mean())
    assert gp["ga"] == pytest.approx(2.06 * np.sqrt(s2g + mse) * gp["h2_bs"] / 100)
    multi = ag.genetic_parameters(_env_data(), ["y"], "trt", replication="block",
                                  environment=["loc", "year"]).table.iloc[0]
    assert 0 <= multi["h2_bs"] <= 100
    assert multi["s2ge"] >= 0


def test_genetic_parameters_negative_component_zero():
    rng = np.random.default_rng(3)
    df = pd.DataFrame({"g": list("ABCD") * 3, "b": np.repeat([1, 2, 3], 4),
                       "y": rng.normal(10, 1, 12)})
    res = ag.genetic_parameters(df, ["y"], "g", replication="b")
    assert res.table["var_g"].iloc[0] >= 0


# ---- correlation / PCA / two sample --------------------------------------

def test_correlation_and_pca(tmp_path):
    rng = np.random.default_rng(0)
    x = rng.normal(size=40)
    df = pd.DataFrame({"a": x, "b": 2 * x + rng.normal(0, .1, 40), "c": rng.normal(size=40)})
    c = ag.correlation(df)
    assert c.r.loc["a", "b"] == pytest.approx(stats.pearsonr(df["a"], df["b"])[0])
    assert c.p.loc["a", "b"] < 0.001
    assert c.stars().loc["a", "b"].endswith("***")
    assert ag.correlation(df, method="spearman").method == "spearman"
    c.plot()
    c.export(tmp_path / "c.xlsx")
    p = ag.pca(df)
    assert p.eigenvalues["variance_pct"].sum() == pytest.approx(100)
    assert p.eigenvalues["eigenvalue"].sum() == pytest.approx(3)
    assert np.allclose((p.loadings**2).sum(axis=1), 1)  # communalities over all PCs
    p.plot()
    with pytest.raises(AnalysisError):
        ag.correlation(df, ["a"])


def test_two_sample_matches_scipy():
    rng = np.random.default_rng(4)
    df = pd.DataFrame({"g": ["a"] * 8 + ["b"] * 10, "y": np.r_[rng.normal(5, 1, 8), rng.normal(6, 2, 10)]})
    a, b = df.y[:8], df.y[8:]
    w = ag.two_sample_test(df, "y", "g")
    ref = stats.ttest_ind(b, a, equal_var=False)
    assert w.p_value == pytest.approx(ref.pvalue)
    assert w.statistic == pytest.approx(ref.statistic)
    assert w.ci[0] < w.estimate < w.ci[1]
    m = ag.two_sample_test(df, "y", "g", test="wilcoxon")
    assert m.p_value == pytest.approx(stats.mannwhitneyu(b, a).pvalue)
    df["id"] = list(range(8)) + [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    with pytest.raises(AnalysisError):
        ag.two_sample_test(df, "y", "g", pair="id")
    p = ag.two_sample_test(df[df.id < 8], "y", "g", test="wilcoxon", pair="id")
    assert "signed-rank" in p.test
    with pytest.raises(AnalysisError):
        ag.two_sample_test(df.assign(g="a"), "y", "g")


# ---- one entry point: ag.analyze ------------------------------------------

def test_analyze_dispatches_alpha_lattice():
    d = _lattice_data()
    res = ag.analyze(d, "y", "alpha_lattice", treatment="treatment", replicate="replicate",
                     block="block", posthoc="dmrt")
    assert isinstance(res, ag.analysis.AlphaLatticeResult)
    assert list(res.comparisons) == ["dmrt"]
    assert "Group" in res.thesis_table()
    none = ag.analyze(d, "y", "alpha_lattice", treatment="treatment", replicate="replicate",
                      block="block", posthoc=None)
    assert none.comparisons == {}
    with pytest.raises(AnalysisError):
        ag.analyze(d, "y", "alpha_lattice", treatment="treatment", block="block")


def test_analyze_dispatches_multi_environment_and_traits(tmp_path):
    df = _env_data()
    res = ag.analyze(df, "y", "rcbd", treatment="trt", block="block",
                     environment=["loc", "year"], posthoc="all")
    assert isinstance(res, ag.analysis.MultiEnvResult)
    assert len(res.comparisons) == 4
    assert sum(c.startswith("Group") for c in res.thesis_table()) == 4
    res.export(tmp_path / "m.xlsx")
    with pytest.raises(AnalysisError):
        ag.analyze(df, "y", "factorial", factors=["trt", "loc"], environment="year")
    both = ag.analyze(df, ["y", "z"], "rcbd", treatment="trt", block="block",
                      environment=["loc", "year"], posthoc="lsd")
    assert isinstance(both, ag.analysis.TraitsReport)
    assert set(both.summary["trait"]) == {"y", "z"}
    assert "bartlett_p" in both.summary
    both.export(tmp_path / "t.xlsx")
    lat = ag.analyze(_lattice_data().assign(w=lambda x: x["y"] * 2), ["y", "w"], "alpha_lattice",
                     treatment="treatment", replicate="replicate", block="block")
    assert list(lat.summary["trait"]) == ["y", "w"]
    plain = ag.analyze(df[(df["loc"] == "A") & (df["year"] == 2023)], ["y", "z"], "rcbd",
                       treatment="trt", block="block")
    assert isinstance(plain, ag.analysis.TraitsReport)


def test_analyze_genetics_flag(tmp_path):
    df = _env_data()
    one = df[(df["loc"] == "A") & (df["year"] == 2023)]
    kw = {"treatment": "trt", "block": "block"}
    rep = ag.analyze(one, "y", "rcbd", genetics=True, **kw)
    direct = ag.genetic_parameters(one, ["y"], "trt", replication="block").table
    assert rep.genetics.table["h2_bs"].iloc[0] == pytest.approx(direct["h2_bs"].iloc[0])
    assert "Genetic parameters" in rep.tables()
    assert "Genetic parameters" in str(rep)
    assert ag.analyze(one, "y", "rcbd", **kw).genetics is None

    multi = ag.analyze(df, "y", "rcbd", genetics=True, environment=["loc", "year"], **kw)
    assert "s2ge" in multi.genetics.table
    multi.export(tmp_path / "g.xlsx")

    many = ag.analyze(one, ["y", "z"], "rcbd", genetics=True, genetic_k=1.76, **kw)
    assert list(many.genetics.table["trait"]) == ["y", "z"]
    assert many.genetics.k == 1.76

    with pytest.raises(AnalysisError):
        ag.analyze(one, "y", "factorial", factors=["trt", "block"], genetics=True)
    with pytest.raises(AnalysisError):
        ag.analyze(_lattice_data(), "y", "alpha_lattice", treatment="treatment",
                   replicate="replicate", block="block", genetics=True)

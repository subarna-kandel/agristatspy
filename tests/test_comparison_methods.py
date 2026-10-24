"""Tests for Tukey (own implementation), LSD, Bonferroni and letters."""
import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import agristatspy as ag
from agristatspy.datasets import load_crd_example, load_rcbd_example
from agristatspy.exceptions import AnalysisError


@pytest.fixture
def crd():
    return load_crd_example()


class TestTukeyAgainstStatsmodels:
    def test_matches_statsmodels_crd(self, crd):
        mine = ag.tukey_hsd(crd, "yield", "treatment").table
        sm = pairwise_tukeyhsd(crd["yield"], crd["treatment"], alpha=0.05)
        ref = pd.DataFrame(sm.summary().data[1:], columns=sm.summary().data[0])
        for col_mine, col_ref in [
            ("mean_diff", "meandiff"), ("p_adj", "p-adj"),
            ("lower_ci", "lower"), ("upper_ci", "upper"),
        ]:
            np.testing.assert_allclose(
                mine[col_mine].astype(float), ref[col_ref].astype(float), atol=2e-3
            )
        assert list(mine["reject"]) == list(ref["reject"])

    def test_unequal_group_sizes(self, crd):
        unbalanced = crd.drop(index=[0, 5])
        mine = ag.tukey_hsd(unbalanced, "yield", "treatment").table
        sm = pairwise_tukeyhsd(unbalanced["yield"], unbalanced["treatment"])
        ref = pd.DataFrame(sm.summary().data[1:], columns=sm.summary().data[0])
        np.testing.assert_allclose(
            mine["p_adj"].astype(float), ref["p-adj"].astype(float), atol=2e-3
        )


class TestBlockedComparisons:
    def test_rcbd_uses_block_error_term(self):
        df = load_rcbd_example()
        crd_res = ag.tukey_hsd(df, "yield", "treatment")
        rcbd_res = ag.tukey_hsd(df, "yield", "treatment", block="block")
        anova = ag.anova_rcbd(df, "yield", "treatment", "block")
        ms_err = anova.table.set_index("source").loc["Error", "ms"]
        assert rcbd_res.mse == pytest.approx(ms_err)
        assert rcbd_res.df_error == 6
        assert rcbd_res.mse < crd_res.mse  # blocking removes block variation

    def test_incomplete_block_layout_rejected(self):
        df = load_rcbd_example().iloc[1:]
        with pytest.raises(AnalysisError):
            ag.tukey_hsd(df, "yield", "treatment", block="block")


class TestLatinSquareComparison:
    def test_uses_latin_square_error_term(self):
        from agristatspy.datasets import load_latin_square_example

        df = load_latin_square_example()
        res = ag.tukey_hsd(df, "yield", "treatment", block=["row", "column"])
        anova = ag.anova_latin_square(df, "yield", "treatment", "row", "column")
        err = anova.table.set_index("source").loc["Error"]
        assert res.mse == pytest.approx(err["ms"])
        assert res.df_error == err["df"] == 6


class TestLSDAndBonferroni:
    def test_lsd_p_equals_unadjusted_t(self, crd):
        res = ag.lsd_test(crd, "yield", "treatment")
        from scipy import stats
        row = res.table.iloc[0]
        g1, g2 = row["group1"], row["group2"]
        d = crd.groupby("treatment")["yield"].mean()
        se = np.sqrt(res.mse * (1 / 4 + 1 / 4))
        p = 2 * stats.t.sf(abs(d[g2] - d[g1]) / se, res.df_error)
        assert row["p_adj"] == pytest.approx(p)

    def test_bonferroni_is_lsd_times_m_capped(self, crd):
        lsd = ag.lsd_test(crd, "yield", "treatment").table["p_adj"].values
        bon = ag.bonferroni(crd, "yield", "treatment").table["p_adj"].values
        np.testing.assert_allclose(bon, np.minimum(1.0, lsd * 3))

    def test_ordering_of_conservativeness(self, crd):
        w = {
            m: (a.upper_ci - a.lower_ci).mean()
            for m, a in {
                "lsd": ag.lsd_test(crd, "yield", "treatment").table,
                "tukey": ag.tukey_hsd(crd, "yield", "treatment").table,
                "bonf": ag.bonferroni(crd, "yield", "treatment").table,
            }.items()
        }
        assert w["lsd"] < w["tukey"] < w["bonf"] * 1.0001

    def test_lsd_warning_in_notes(self, crd):
        assert "NOT controlled" in ag.lsd_test(crd, "yield", "treatment").model_notes

    @pytest.mark.parametrize("alpha", [0, 1, -0.1, 2, True, "0.05"])
    def test_invalid_alpha(self, crd, alpha):
        with pytest.raises(AnalysisError):
            ag.tukey_hsd(crd, "yield", "treatment", alpha=alpha)


class TestLetters:
    def test_clear_separation_gives_distinct_letters(self, crd):
        letters = ag.tukey_hsd(crd, "yield", "treatment").letters()
        assert list(letters["treatment"]) == ["N100", "N50", "Control"]
        assert list(letters["group"]) == ["a", "b", "c"]

    def test_no_difference_shares_one_letter(self):
        rng = np.random.default_rng(0)
        df = pd.DataFrame(
            {"t": list("ABC") * 5, "y": rng.normal(10, 1, 15)}
        )
        letters = ag.tukey_hsd(df, "y", "t").letters()
        assert set(letters["group"]) == {"a"}

    def test_overlapping_groups(self):
        # A and C differ; B overlaps both -> a, ab, b
        noise = [-0.8, -0.3, 0.3, 0.8]
        df = pd.DataFrame(
            {
                "t": ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
                "y": [10 + e for e in noise] + [11 + e for e in noise] + [12 + e for e in noise],
            }
        )
        res = ag.lsd_test(df, "y", "t", alpha=0.05)
        tbl = res.table.set_index(["group1", "group2"])["reject"]
        assert bool(tbl[("A", "C")])
        assert not bool(tbl[("A", "B")])
        assert not bool(tbl[("B", "C")])
        g = dict(zip(*[res.letters()[c] for c in ("treatment", "group")], strict=True))
        assert g == {"C": "a", "B": "ab", "A": "b"}

    def test_letters_consistent_with_pairwise(self, crd):
        res = ag.lsd_test(crd, "yield", "treatment")
        g = dict(zip(res.letters()["treatment"], res.letters()["group"], strict=True))
        for _, r in res.table.iterrows():
            share = bool(set(g[r["group1"]]) & set(g[r["group2"]]))
            assert share == (not r["reject"])

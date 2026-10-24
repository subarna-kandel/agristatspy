"""Tests for factorial designs and ANOVA."""
import pandas as pd
import pytest

import agristatspy as ag
from agristatspy.datasets import load_factorial_example
from agristatspy.exceptions import AnalysisError, DesignError

FACTORS = {"N": ["N0", "N50", "N100"], "V": ["V1", "V2"]}


class TestFactorialDesign:
    def test_crd_layout_counts(self):
        d = ag.design.factorial(FACTORS, replications=3, layout="crd", seed=1)
        assert len(d.table) == 18
        assert (d.table.groupby(["N", "V"]).size() == 3).all()
        assert d.metadata.design_type == "Factorial-CRD"

    def test_rcbd_layout_each_combo_once_per_block(self):
        d = ag.design.factorial(FACTORS, replications=4, layout="rcbd", seed=1)
        assert len(d.table) == 24
        assert (d.table.groupby(["block", "N", "V"]).size() == 1).all()

    def test_seed_reproducible(self):
        a = ag.design.factorial(FACTORS, replications=2, seed=5).table
        b = ag.design.factorial(FACTORS, replications=2, seed=5).table
        assert a.equals(b)

    def test_text_layout_and_header(self):
        d = ag.design.factorial(FACTORS, replications=2, layout="rcbd", seed=1)
        assert d.text_layout().startswith("Factorial-RCBD")
        assert "BLOCK 2" in d.text_layout()

    @pytest.mark.parametrize(
        "bad",
        [{}, {"A": ["a1", "a2"]}, {"A": ["a1"], "B": ["b1", "b2"]}, {"A": ["x", "x"], "B": ["b1", "b2"]}, "AB"],
    )
    def test_invalid_factors(self, bad):
        with pytest.raises(DesignError):
            ag.design.factorial(bad, replications=2)

    def test_invalid_layout(self):
        with pytest.raises(DesignError):
            ag.design.factorial(FACTORS, layout="latin")


class TestFactorialANOVA:
    @pytest.fixture
    def df(self):
        return load_factorial_example()

    def test_rcbd_table_structure(self, df):
        r = ag.anova_factorial(df, "yield", ["variety", "nitrogen"], block="block")
        dfs = dict(zip(r.table["source"], r.table["df"], strict=True))
        assert dfs == {
            "Block": 2, "variety": 1, "nitrogen": 2, "variety:nitrogen": 2,
            "Error": 10, "Total": 17,
        }

    def test_ss_partition(self, df):
        r = ag.anova_factorial(df, "yield", ["variety", "nitrogen"], block="block")
        t = r.table.set_index("source")["ss"]
        assert t.drop("Total").sum() == pytest.approx(t["Total"])

    def test_hand_check_main_effect_ss(self, df):
        r = ag.anova_factorial(df, "yield", ["variety", "nitrogen"])
        g = df["yield"].mean()
        expected = sum(
            len(x) * (x.mean() - g) ** 2 for _, x in df.groupby("nitrogen")["yield"]
        )
        assert r.table.set_index("source").loc["nitrogen", "ss"] == pytest.approx(expected)

    def test_crd_dispatch_matches(self, df):
        a = ag.anova_factorial(df, "yield", ["variety", "nitrogen"])
        b = ag.anova(df, "yield", "factorial", factors=["variety", "nitrogen"])
        pd.testing.assert_frame_equal(a.table, b.table)

    def test_unbalanced_rejected(self, df):
        with pytest.raises(AnalysisError):
            ag.anova_factorial(df.iloc[1:], "yield", ["variety", "nitrogen"])

    def test_missing_response_rejected_not_silently_used(self, df):
        d = df.copy()
        d.loc[0, "yield"] = float("nan")
        with pytest.raises(AnalysisError):
            ag.anova_factorial(d, "yield", ["variety", "nitrogen"], block="block")

    def test_single_factor_rejected(self, df):
        with pytest.raises(AnalysisError):
            ag.anova_factorial(df, "yield", ["variety"])

    def test_no_error_df_rejected(self):
        d = ag.design.factorial(FACTORS, replications=1, seed=1).table
        d["y"] = range(len(d))
        with pytest.raises(AnalysisError):
            ag.anova_factorial(d, "y", ["N", "V"])

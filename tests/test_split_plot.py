"""Tests for the split-plot design and ANOVA."""
import pytest

import agristatspy as ag
from agristatspy.datasets import load_split_plot_example
from agristatspy.exceptions import AnalysisError, DesignError


class TestSplitPlotDesign:
    @pytest.fixture
    def d(self):
        return ag.design.split_plot(
            ["I0", "I1", "I2"], ["V1", "V2", "V3", "V4"], replications=3,
            main_name="irr", sub_name="var", seed=7,
        )

    def test_counts(self, d):
        assert len(d.table) == 36
        assert d.metadata.total_plots == 36

    def test_each_main_level_once_per_block(self, d):
        t = d.table
        assert (t.groupby(["block", "irr"])["main_plot"].nunique() == 1).all()
        assert (t.groupby("block")["irr"].nunique() == 3).all()

    def test_each_sub_level_once_per_main_plot(self, d):
        t = d.table
        assert (t.groupby(["main_plot", "var"]).size() == 1).all()
        assert (t.groupby("main_plot").size() == 4).all()

    def test_plot_ids_unique(self, d):
        assert d.table["plot_id"].is_unique

    def test_seed_reproducible(self, d):
        again = ag.design.split_plot(
            ["I0", "I1", "I2"], ["V1", "V2", "V3", "V4"], replications=3,
            main_name="irr", sub_name="var", seed=7,
        )
        assert d.table.equals(again.table)

    def test_text_layout(self, d):
        text = d.text_layout()
        assert "Split-plot" in text
        assert "BLOCK 3" in text
        assert "Main plot B1-M1" in text
        for pid in d.table["plot_id"]:
            assert pid in text

    @pytest.mark.parametrize(
        "args",
        [
            (["A"], ["x", "y"]),
            (["A", "A"], ["x", "y"]),
            (["A", "B"], "xy"),
        ],
    )
    def test_invalid_levels(self, args):
        with pytest.raises(DesignError):
            ag.design.split_plot(*args, replications=2)

    def test_name_clash(self):
        with pytest.raises(DesignError):
            ag.design.split_plot(["a", "b"], ["x", "y"], main_name="block", sub_name="s")


class TestSplitPlotANOVA:
    @pytest.fixture
    def result(self):
        return ag.anova_split_plot(
            load_split_plot_example(), "yield", "irrigation", "variety", "block"
        )

    def test_degrees_of_freedom(self, result):
        dfs = dict(zip(result.table["source"], result.table["df"], strict=True))
        assert dfs == {
            "Block": 2, "irrigation": 1, "Error (a)": 2, "variety": 2,
            "irrigation:variety": 2, "Error (b)": 8, "Total": 17,
        }

    def test_ss_partition(self, result):
        t = result.table.set_index("source")["ss"]
        assert t.drop("Total").sum() == pytest.approx(t["Total"])

    def test_error_a_hand_check(self):
        df = load_split_plot_example()
        g = df["yield"].mean()
        ss_blk = sum(len(x) * (x.mean() - g) ** 2 for _, x in df.groupby("block")["yield"])
        ss_irr = sum(len(x) * (x.mean() - g) ** 2 for _, x in df.groupby("irrigation")["yield"])
        cell = df.groupby(["block", "irrigation"])["yield"].mean()
        ss_whole = 3 * ((cell - g) ** 2).sum()
        expected = ss_whole - ss_blk - ss_irr
        r = ag.anova_split_plot(df, "yield", "irrigation", "variety", "block")
        assert r.table.set_index("source").loc["Error (a)", "ss"] == pytest.approx(expected)

    def test_main_tested_against_error_a(self, result):
        t = result.table.set_index("source")
        assert t.loc["irrigation", "f_stat"] == pytest.approx(
            t.loc["irrigation", "ms"] / t.loc["Error (a)", "ms"]
        )
        assert t.loc["variety", "f_stat"] == pytest.approx(
            t.loc["variety", "ms"] / t.loc["Error (b)", "ms"]
        )

    def test_dispatcher(self):
        r = ag.anova(
            load_split_plot_example(), "yield", "split_plot",
            main="irrigation", sub="variety", block="block",
        )
        assert r.design == "Split-plot"

    def test_incomplete_rejected(self):
        df = load_split_plot_example().iloc[1:]
        with pytest.raises(AnalysisError):
            ag.anova_split_plot(df, "yield", "irrigation", "variety", "block")

    def test_dispatcher_unknown_and_missing_args(self):
        df = load_split_plot_example()
        with pytest.raises(AnalysisError):
            ag.anova(df, "yield", "mystery")
        with pytest.raises(AnalysisError):
            ag.anova(df, "yield", "split_plot", main="irrigation")

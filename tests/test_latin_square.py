"""Tests for the Latin square design and its ANOVA."""
import pytest

import agristatspy as ag
from agristatspy.datasets import load_latin_square_example
from agristatspy.exceptions import AnalysisError, DesignError


class TestLatinSquareDesign:
    @pytest.mark.parametrize("t", [2, 3, 4, 5, 8])
    def test_each_treatment_once_per_row_and_column(self, t):
        d = ag.design.latin_square(n_treatments=t, seed=t)
        tbl = d.table
        assert len(tbl) == t * t
        for axis in ("row", "column"):
            assert (tbl.groupby([axis, "treatment"]).size() == 1).all()

    def test_seed_reproducible(self):
        a = ag.design.latin_square(n_treatments=5, seed=9).table
        b = ag.design.latin_square(n_treatments=5, seed=9).table
        assert a.equals(b)

    def test_different_seeds_vary(self):
        tables = {
            tuple(ag.design.latin_square(n_treatments=5, seed=s).table["treatment"])
            for s in range(10)
        }
        assert len(tables) > 1

    def test_labels_and_validate(self):
        d = ag.design.latin_square(treatments=["a", "b", "c"], seed=1)
        assert set(d.table["treatment"]) == {"a", "b", "c"}
        assert d.validate() is True

    def test_single_treatment_rejected(self):
        with pytest.raises(DesignError):
            ag.design.latin_square(n_treatments=1)

    def test_text_and_plot_layout(self):
        d = ag.design.latin_square(n_treatments=3, seed=1)
        text = d.text_layout()
        assert "Latin square" in text
        assert "Row 3" in text
        assert "Col 3" in text
        for pid in d.table["plot_id"]:
            assert pid in text


class TestLatinSquareANOVA:
    @pytest.fixture
    def result(self):
        df = load_latin_square_example()
        return ag.anova_latin_square(df, "yield", "treatment", "row", "column")

    def test_degrees_of_freedom(self, result):
        dfs = dict(zip(result.table["source"], result.table["df"], strict=True))
        assert dfs == {"Treatment": 3, "Row": 3, "Column": 3, "Error": 6, "Total": 15}

    def test_ss_partition(self, result):
        t = result.table.set_index("source")["ss"]
        assert t["Treatment"] + t["Row"] + t["Column"] + t["Error"] == pytest.approx(t["Total"])

    def test_treatment_effect_detected(self, result):
        p = result.table.set_index("source").loc["Treatment", "p_value"]
        assert p < 0.001

    def test_dispatcher(self):
        df = load_latin_square_example()
        r = ag.anova(df, "yield", "latin_square", treatment="treatment", row="row", column="column")
        assert r.design == "Latin square"

    def test_nonsquare_rejected(self):
        df = load_latin_square_example()
        with pytest.raises(AnalysisError):
            ag.anova_latin_square(df[df["row"] != 4], "yield", "treatment", "row", "column")

    def test_non_latin_rejected(self):
        df = load_latin_square_example().copy()
        df.loc[0, "treatment"] = "N50"  # breaks the once-per-row/column rule
        with pytest.raises(AnalysisError):
            ag.anova_latin_square(df, "yield", "treatment", "row", "column")

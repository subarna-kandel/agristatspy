"""
Tests for Tukey HSD mean comparisons.
"""
import pandas as pd
import pytest

import agristatspy as ag
from agristatspy.exceptions import AnalysisError

DATA = pd.DataFrame(
    {
        "treatment": ["Control"] * 4 + ["N50"] * 4 + ["N100"] * 4,
        "yield": [2.3, 2.7, 2.5, 2.4, 3.4, 3.6, 3.5, 3.3, 4.4, 4.6, 4.5, 4.5],
    }
)


class TestTukeyHSD:
    @pytest.fixture
    def result(self):
        return ag.tukey_hsd(DATA, response="yield", treatment="treatment")

    def test_method_label(self, result):
        assert result.method == "Tukey HSD"

    def test_alpha_stored(self, result):
        assert result.alpha == 0.05

    def test_table_columns(self, result):
        expected = {"group1", "group2", "mean_diff", "lower_ci", "upper_ci", "p_adj", "reject"}
        assert set(result.table.columns) == expected

    def test_number_of_comparisons(self, result):
        # C(3,2) = 3 pairwise comparisons
        assert len(result.table) == 3

    def test_control_vs_n100_significant(self, result):
        """Control vs N100 difference is large and should be significant."""
        tbl = result.table
        pair = tbl[
            ((tbl["group1"] == "Control") & (tbl["group2"] == "N100")) |
            ((tbl["group1"] == "N100") & (tbl["group2"] == "Control"))
        ]
        assert len(pair) == 1
        assert bool(pair["reject"].iloc[0]) is True

    def test_str_representation(self, result):
        s = str(result)
        assert "Tukey" in s

    def test_missing_column_raises(self):
        with pytest.raises(AnalysisError):
            ag.tukey_hsd(DATA, response="nonexistent", treatment="treatment")

    def test_single_group_raises(self):
        single = DATA[DATA["treatment"] == "Control"]
        with pytest.raises(AnalysisError):
            ag.tukey_hsd(single, response="yield", treatment="treatment")

"""
Tests for descriptive statistics (treatment_summary).
"""
import math

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


class TestTreatmentSummary:
    @pytest.fixture
    def summary(self):
        return ag.treatment_summary(DATA, response="yield", treatment="treatment")

    def test_three_rows(self, summary):
        assert len(summary) == 3

    def test_columns_present(self, summary):
        for col in ["treatment", "n", "mean", "std", "se", "min", "median", "max", "cv_pct"]:
            assert col in summary.columns, f"Column {col!r} missing"

    def test_n_per_treatment(self, summary):
        for _, row in summary.iterrows():
            assert row["n"] == 4

    def test_control_mean(self, summary):
        control = summary[summary["treatment"] == "Control"].iloc[0]
        assert abs(control["mean"] - 2.475) < 1e-9

    def test_n100_mean(self, summary):
        n100 = summary[summary["treatment"] == "N100"].iloc[0]
        assert abs(n100["mean"] - 4.5) < 1e-9

    def test_cv_pct_positive(self, summary):
        cv = summary["cv_pct"].iloc[0]
        assert cv > 0

    def test_cv_note_in_attrs(self, summary):
        assert "cv_note" in summary.attrs

    def test_missing_response_raises(self):
        with pytest.raises(AnalysisError):
            ag.treatment_summary(DATA, response="no_col", treatment="treatment")

    def test_non_numeric_response_raises(self):
        bad = DATA.copy()
        bad["yield"] = "text"
        with pytest.raises(AnalysisError):
            ag.treatment_summary(bad, response="yield", treatment="treatment")

    def test_zero_mean_cv_is_nan(self):
        """When grand mean is zero, cv_pct should be NaN."""
        zero_data = pd.DataFrame(
            {
                "treatment": ["A"] * 3 + ["B"] * 3,
                "yield": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            }
        )
        result = ag.treatment_summary(zero_data, response="yield", treatment="treatment")
        assert math.isnan(result["cv_pct"].iloc[0])

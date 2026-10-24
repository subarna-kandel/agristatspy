"""
Tests for ANOVA (CRD and RCBD).

All expected values are hand-calculated from first principles or verified
against scipy.stats.f_oneway and statsmodels.  Floating-point assertions use
tolerances appropriate for double-precision arithmetic.
"""

import numpy as np
import pandas as pd
import pytest
import scipy.stats as sp_stats

import agristatspy as ag
from agristatspy.exceptions import AnalysisError

# ---------------------------------------------------------------------------
# Shared fixture data
# ---------------------------------------------------------------------------

# CRD: 3 treatments × 4 reps — yields that produce a clearly significant F.
CRD_DATA = pd.DataFrame(
    {
        "treatment": ["Control"] * 4 + ["N50"] * 4 + ["N100"] * 4,
        "yield": [2.3, 2.7, 2.5, 2.4, 3.4, 3.6, 3.5, 3.3, 4.4, 4.6, 4.5, 4.5],
    }
)

# RCBD: 4 treatments × 3 blocks.
RCBD_DATA = pd.DataFrame(
    {
        "treatment": ["Control", "N50", "N100", "N150"] * 3,
        "block": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        "yield": [2.1, 3.2, 4.1, 4.8, 2.3, 3.5, 4.3, 5.0, 2.0, 3.1, 4.0, 4.7],
    }
)


# ---------------------------------------------------------------------------
# CRD ANOVA
# ---------------------------------------------------------------------------

class TestANOVACRD:
    @pytest.fixture
    def result(self):
        return ag.anova_crd(CRD_DATA, response="yield", treatment="treatment")

    def test_design_label(self, result):
        assert result.design == "CRD"

    def test_n_obs(self, result):
        assert result.n_obs == 12

    def test_grand_mean(self, result):
        assert abs(result.grand_mean - np.mean(CRD_DATA["yield"])) < 1e-9

    def test_df_treatment(self, result):
        trt_row = result.table[result.table["source"] == "Treatment"].iloc[0]
        assert int(trt_row["df"]) == 2  # k-1 = 3-1

    def test_df_error(self, result):
        err_row = result.table[result.table["source"] == "Error"].iloc[0]
        assert int(err_row["df"]) == 9  # n-k = 12-3

    def test_df_total(self, result):
        tot_row = result.table[result.table["source"] == "Total"].iloc[0]
        assert int(tot_row["df"]) == 11  # n-1

    def test_ss_partition(self, result):
        """SS_treatment + SS_error must equal SS_total."""
        trt_ss = float(result.table[result.table["source"] == "Treatment"]["ss"].iloc[0])
        err_ss = float(result.table[result.table["source"] == "Error"]["ss"].iloc[0])
        tot_ss = float(result.table[result.table["source"] == "Total"]["ss"].iloc[0])
        assert abs(trt_ss + err_ss - tot_ss) < 1e-9

    def test_f_stat_matches_scipy(self, result):
        """F-statistic must match scipy.stats.f_oneway."""
        groups = [
            CRD_DATA[CRD_DATA["treatment"] == t]["yield"].values
            for t in ["Control", "N50", "N100"]
        ]
        f_scipy, _ = sp_stats.f_oneway(*groups)
        trt_row = result.table[result.table["source"] == "Treatment"].iloc[0]
        assert abs(float(trt_row["f_stat"]) - f_scipy) < 1e-6

    def test_p_value_significant(self, result):
        """Treatment effect in this dataset should be highly significant."""
        trt_row = result.table[result.table["source"] == "Treatment"].iloc[0]
        assert float(trt_row["p_value"]) < 0.001

    def test_cv_pct_positive(self, result):
        assert result.cv_pct > 0

    def test_str_representation(self, result):
        s = str(result)
        assert "ANOVA" in s
        assert "CRD" in s

    def test_missing_response_column_raises(self):
        with pytest.raises(AnalysisError):
            ag.anova_crd(CRD_DATA, response="nonexistent", treatment="treatment")

    def test_missing_treatment_column_raises(self):
        with pytest.raises(AnalysisError):
            ag.anova_crd(CRD_DATA, response="yield", treatment="nonexistent")


# ---------------------------------------------------------------------------
# RCBD ANOVA
# ---------------------------------------------------------------------------

class TestANOVARCBD:
    @pytest.fixture
    def result(self):
        return ag.anova_rcbd(RCBD_DATA, response="yield", treatment="treatment", block="block")

    def test_design_label(self, result):
        assert result.design == "RCBD"

    def test_n_obs(self, result):
        assert result.n_obs == 12

    def test_df_treatment(self, result):
        trt_row = result.table[result.table["source"] == "Treatment"].iloc[0]
        assert int(trt_row["df"]) == 3  # k-1 = 4-1

    def test_df_block(self, result):
        blk_row = result.table[result.table["source"] == "Block"].iloc[0]
        assert int(blk_row["df"]) == 2  # b-1 = 3-1

    def test_df_error(self, result):
        err_row = result.table[result.table["source"] == "Error"].iloc[0]
        assert int(err_row["df"]) == 6  # (k-1)(b-1) = 3×2

    def test_df_total(self, result):
        tot_row = result.table[result.table["source"] == "Total"].iloc[0]
        assert int(tot_row["df"]) == 11

    def test_ss_partition(self, result):
        trt_ss = float(result.table[result.table["source"] == "Treatment"]["ss"].iloc[0])
        blk_ss = float(result.table[result.table["source"] == "Block"]["ss"].iloc[0])
        err_ss = float(result.table[result.table["source"] == "Error"]["ss"].iloc[0])
        tot_ss = float(result.table[result.table["source"] == "Total"]["ss"].iloc[0])
        assert abs(trt_ss + blk_ss + err_ss - tot_ss) < 1e-6

    def test_treatment_f_significant(self, result):
        trt_row = result.table[result.table["source"] == "Treatment"].iloc[0]
        assert float(trt_row["p_value"]) < 0.001

    def test_cv_pct_positive(self, result):
        assert result.cv_pct > 0

    def test_missing_block_column_raises(self):
        with pytest.raises(AnalysisError):
            ag.anova_rcbd(RCBD_DATA, response="yield", treatment="treatment", block="no_block")

    def test_incomplete_blocks_raises(self):
        """Data with duplicate treatment-block combinations should raise AnalysisError."""
        bad = pd.concat([RCBD_DATA, RCBD_DATA.iloc[[0]]], ignore_index=True)
        with pytest.raises(AnalysisError):
            ag.anova_rcbd(bad, response="yield", treatment="treatment", block="block")

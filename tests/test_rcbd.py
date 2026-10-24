"""
Tests for RCBD (Randomized Complete Block Design) generation.

Covers:
- Total plots (k × b)
- Number of blocks
- Complete-block invariant (every treatment in every block, exactly once)
- Treatment frequencies
- Plot-ID uniqueness
- Position within block values
- Seed reproducibility
- Validation errors
"""
import pytest

import agristatspy as ag
from agristatspy.exceptions import DesignError

# ---------------------------------------------------------------------------
# Basic invariants: 4 treatments × 3 blocks = 12 plots
# ---------------------------------------------------------------------------

class TestRCBDInvariants:
    @pytest.fixture
    def design_4x3(self):
        return ag.design.rcbd(
            treatments=["Control", "N50", "N100", "N150"],
            replications=3,
            seed=42,
        )

    def test_total_plots(self, design_4x3):
        assert len(design_4x3.table) == 12

    def test_number_of_blocks(self, design_4x3):
        assert design_4x3.table["block"].nunique() == 3

    def test_plot_id_uniqueness(self, design_4x3):
        assert design_4x3.table["plot_id"].nunique() == 12

    def test_every_treatment_in_every_block(self, design_4x3):
        tbl = design_4x3.table
        for block_id in tbl["block"].unique():
            block_treatments = set(tbl[tbl["block"] == block_id]["treatment"].tolist())
            assert block_treatments == {"Control", "N50", "N100", "N150"}, (
                f"Block {block_id} does not contain all treatments: {block_treatments}"
            )

    def test_each_treatment_once_per_block(self, design_4x3):
        tbl = design_4x3.table
        for block_id in tbl["block"].unique():
            block_trt = tbl[tbl["block"] == block_id]["treatment"]
            assert len(block_trt) == 4
            assert block_trt.nunique() == 4

    def test_each_treatment_appears_n_blocks_times(self, design_4x3):
        counts = design_4x3.table.groupby("treatment").size()
        for trt in ["Control", "N50", "N100", "N150"]:
            assert counts[trt] == 3

    def test_position_in_block_values(self, design_4x3):
        tbl = design_4x3.table
        for block_id in tbl["block"].unique():
            positions = sorted(tbl[tbl["block"] == block_id]["position_in_block"].tolist())
            assert positions == [1, 2, 3, 4]


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

class TestRCBDSeedReproducibility:
    def test_same_seed_same_allocation(self):
        d1 = ag.design.rcbd(n_treatments=4, replications=3, seed=42)
        d2 = ag.design.rcbd(n_treatments=4, replications=3, seed=42)
        assert d1.table["treatment"].tolist() == d2.table["treatment"].tolist()

    def test_different_seeds_different_allocations(self):
        d1 = ag.design.rcbd(n_treatments=4, replications=3, seed=1)
        d2 = ag.design.rcbd(n_treatments=4, replications=3, seed=2)
        assert d1.table["treatment"].tolist() != d2.table["treatment"].tolist()

    def test_none_seed_runs_without_error(self):
        d = ag.design.rcbd(n_treatments=4, replications=3, seed=None)
        assert len(d.table) == 12


# ---------------------------------------------------------------------------
# Auto-generated labels
# ---------------------------------------------------------------------------

class TestRCBDAutoLabels:
    def test_auto_labels_format(self):
        d = ag.design.rcbd(n_treatments=3, replications=2, seed=0)
        assert set(d.table["treatment"].unique()) == {"T1", "T2", "T3"}


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

class TestRCBDMetadata:
    def test_design_type(self):
        d = ag.design.rcbd(n_treatments=3, replications=2, seed=0)
        assert d.metadata.design_type == "RCBD"

    def test_metadata_total_plots(self):
        d = ag.design.rcbd(n_treatments=4, replications=3, seed=0)
        assert d.metadata.total_plots == 12

    def test_metadata_n_blocks(self):
        d = ag.design.rcbd(n_treatments=4, replications=3, seed=0)
        assert d.metadata.extra["n_blocks"] == 3


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

class TestRCBDValidate:
    def test_validate_passes(self):
        d = ag.design.rcbd(n_treatments=4, replications=3, seed=42)
        assert d.validate() is True


# ---------------------------------------------------------------------------
# Input errors
# ---------------------------------------------------------------------------

class TestRCBDInputErrors:
    def test_no_treatments_or_n_treatments_raises(self):
        with pytest.raises(DesignError):
            ag.design.rcbd(replications=3)

    def test_contradictory_n_treatments_raises(self):
        with pytest.raises(DesignError):
            ag.design.rcbd(treatments=["A", "B", "C"], n_treatments=4, replications=3)

    def test_duplicate_treatment_raises(self):
        with pytest.raises(DesignError):
            ag.design.rcbd(treatments=["A", "A", "B"], replications=2)

    def test_zero_replications_raises(self):
        with pytest.raises(DesignError):
            ag.design.rcbd(n_treatments=3, replications=0)

    def test_bool_replications_raises(self):
        with pytest.raises(DesignError):
            ag.design.rcbd(n_treatments=3, replications=True)

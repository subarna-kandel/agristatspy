"""
Tests for CRD (Completely Randomized Design) generation.

Covers:
- Total experimental units (k × r)
- Treatment frequency invariants
- Plot-ID uniqueness
- Seed reproducibility
- Different seed produces different allocation
- Auto-generated treatment labels
- Validation errors for bad inputs
"""
import pytest

import agristatspy as ag
from agristatspy.exceptions import DesignError

# ---------------------------------------------------------------------------
# Basic invariants
# ---------------------------------------------------------------------------

class TestCRDInvariants:
    def test_total_plots(self):
        d = ag.design.crd(treatments=["A", "B", "C", "D"], replications=5, seed=1)
        assert len(d.table) == 20

    def test_treatment_frequencies(self):
        d = ag.design.crd(treatments=["A", "B", "C", "D"], replications=5, seed=1)
        counts = d.table.groupby("treatment").size()
        for t in ["A", "B", "C", "D"]:
            assert counts[t] == 5, f"Treatment {t} expected 5 reps, got {counts[t]}"

    def test_plot_id_uniqueness(self):
        d = ag.design.crd(n_treatments=4, replications=5, seed=1)
        assert d.table["plot_id"].nunique() == 20

    def test_replication_numbers_per_treatment(self):
        """Each treatment should have replication_number 1..r."""
        d = ag.design.crd(treatments=["A", "B", "C"], replications=4, seed=99)
        for trt in ["A", "B", "C"]:
            reps = sorted(d.table[d.table["treatment"] == trt]["replication_number"].tolist())
            assert reps == [1, 2, 3, 4], f"Unexpected rep numbers for {trt}: {reps}"


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

class TestCRDSeedReproducibility:
    def test_same_seed_same_allocation(self):
        d1 = ag.design.crd(n_treatments=4, replications=5, seed=42)
        d2 = ag.design.crd(n_treatments=4, replications=5, seed=42)
        assert d1.table["treatment"].tolist() == d2.table["treatment"].tolist()

    def test_different_seeds_different_allocations(self):
        """Different seeds should (almost always) produce different orderings."""
        d1 = ag.design.crd(n_treatments=4, replications=5, seed=1)
        d2 = ag.design.crd(n_treatments=4, replications=5, seed=2)
        # With k=4, r=5, the probability of identical order is negligible.
        assert d1.table["treatment"].tolist() != d2.table["treatment"].tolist()

    def test_none_seed_runs_without_error(self):
        d = ag.design.crd(n_treatments=3, replications=3, seed=None)
        assert len(d.table) == 9


# ---------------------------------------------------------------------------
# Auto-generated labels
# ---------------------------------------------------------------------------

class TestCRDAutoLabels:
    def test_auto_labels_format(self):
        d = ag.design.crd(n_treatments=4, replications=1, seed=0)
        assert set(d.table["treatment"].unique()) == {"T1", "T2", "T3", "T4"}

    def test_explicit_labels_used(self):
        d = ag.design.crd(treatments=["Control", "Low", "Medium", "High"], replications=2, seed=0)
        assert set(d.table["treatment"].unique()) == {"Control", "Low", "Medium", "High"}


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

class TestCRDMetadata:
    def test_design_type(self):
        d = ag.design.crd(n_treatments=3, replications=3, seed=0)
        assert d.metadata.design_type == "CRD"

    def test_metadata_total_plots(self):
        d = ag.design.crd(n_treatments=3, replications=3, seed=0)
        assert d.metadata.total_plots == 9

    def test_metadata_seed(self):
        d = ag.design.crd(n_treatments=3, replications=3, seed=77)
        assert d.metadata.seed == 77


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestCRDValidate:
    def test_validate_passes(self):
        d = ag.design.crd(n_treatments=4, replications=5, seed=1)
        assert d.validate() is True


# ---------------------------------------------------------------------------
# Input error handling
# ---------------------------------------------------------------------------

class TestCRDInputErrors:
    def test_no_treatments_or_n_treatments_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(replications=3)

    def test_contradictory_n_treatments_and_labels_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(treatments=["A", "B", "C"], n_treatments=2, replications=3)

    def test_duplicate_treatment_labels_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(treatments=["A", "A", "B"], replications=3)

    def test_empty_label_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(treatments=["A", "", "B"], replications=3)

    def test_zero_replications_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=3, replications=0)

    def test_negative_replications_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=3, replications=-1)

    def test_bool_replications_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=3, replications=True)

    def test_zero_n_treatments_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=0, replications=3)

    def test_bool_seed_raises(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=3, replications=2, seed=True)


# ---------------------------------------------------------------------------
# Plot prefix / start_plot
# ---------------------------------------------------------------------------

class TestCRDPlotIDs:
    def test_custom_prefix(self):
        d = ag.design.crd(n_treatments=2, replications=2, seed=0, plot_prefix="PLOT")
        assert all(pid.startswith("PLOT") for pid in d.table["plot_id"])

    def test_start_plot_offset(self):
        d = ag.design.crd(n_treatments=2, replications=2, seed=0, start_plot=100)
        first_num = int(d.table["plot_id"].iloc[0].lstrip("P"))
        assert first_num == 100

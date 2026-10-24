"""
Tests for text layout rendering (CRD and RCBD).

Verifies that:
- The rendered string contains the correct header information.
- All plot IDs and treatment labels from the table appear in the layout.
- RCBD layouts contain block headers and separators.
- Layout rendering does not alter the underlying assignment table.
- Return type is str.
"""
import agristatspy as ag


class TestCRDTextLayout:
    def setup_method(self):
        self.design = ag.design.crd(
            treatments=["Control", "N50", "N100"],
            replications=4,
            seed=12,
        )

    def test_returns_str(self):
        layout = self.design.text_layout()
        assert isinstance(layout, str)

    def test_header_contains_design_type(self):
        layout = self.design.text_layout()
        assert "CRD" in layout

    def test_header_contains_counts(self):
        layout = self.design.text_layout()
        assert "3" in layout   # treatments
        assert "4" in layout   # replications
        assert "12" in layout  # total plots

    def test_all_plot_ids_in_layout(self):
        layout = self.design.text_layout()
        for pid in self.design.table["plot_id"]:
            assert pid in layout, f"Plot ID {pid!r} missing from layout"

    def test_all_treatments_in_layout(self):
        layout = self.design.text_layout()
        for trt in ["Control", "N50", "N100"]:
            assert trt in layout, f"Treatment {trt!r} missing from layout"

    def test_layout_does_not_mutate_table(self):
        before = self.design.table.copy()
        _ = self.design.text_layout()
        after = self.design.table
        assert before.equals(after)

    def test_str_dunder_returns_layout(self):
        assert str(self.design) == self.design.text_layout()


class TestRCBDTextLayout:
    def setup_method(self):
        self.design = ag.design.rcbd(
            treatments=["Control", "N50", "N100", "N150"],
            replications=3,
            seed=42,
        )

    def test_returns_str(self):
        layout = self.design.text_layout()
        assert isinstance(layout, str)

    def test_header_contains_rcbd(self):
        layout = self.design.text_layout()
        assert "RCBD" in layout

    def test_header_contains_counts(self):
        layout = self.design.text_layout()
        assert "4" in layout   # treatments
        assert "3" in layout   # blocks
        assert "12" in layout  # total plots

    def test_all_block_headers_present(self):
        layout = self.design.text_layout()
        for b in [1, 2, 3]:
            assert f"BLOCK {b}" in layout, f"Block {b} header missing"

    def test_block_separator_present_by_default(self):
        layout = self.design.text_layout(block_separator=True, separator="=")
        # At least one line of '=' characters should appear.
        assert any(line.startswith("=") for line in layout.splitlines())

    def test_block_separator_absent_when_disabled(self):
        layout = self.design.text_layout(block_separator=False)
        # No separator-only lines (lines that are all '=' chars).
        separator_lines = [
            line for line in layout.splitlines()
            if line and all(c == "=" for c in line.strip())
        ]
        assert len(separator_lines) == 0

    def test_all_plot_ids_in_layout(self):
        layout = self.design.text_layout()
        for pid in self.design.table["plot_id"]:
            assert pid in layout

    def test_all_treatments_in_layout(self):
        layout = self.design.text_layout()
        for trt in ["Control", "N50", "N100", "N150"]:
            assert trt in layout

    def test_layout_matches_table(self):
        """Every treatment-plot pair from the table must appear in the layout."""
        layout = self.design.text_layout()
        tbl = self.design.table
        for _, row in tbl.iterrows():
            assert row["plot_id"] in layout
            assert row["treatment"] in layout

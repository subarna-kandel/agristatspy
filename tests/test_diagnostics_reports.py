"""Tests for diagnostics, plots, reports and the extra design options."""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import agristatspy as ag  # noqa: E402
from agristatspy.datasets import load_crd_example, load_rcbd_example  # noqa: E402
from agristatspy.exceptions import AnalysisError, DesignError, ExportError  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


class TestDiagnostics:
    def test_crd_residuals_sum_to_zero_within_groups(self):
        d = ag.residual_diagnostics(load_crd_example(), "yield", "treatment")
        sums = d.residuals.groupby("treatment")["residual"].sum()
        assert sums.abs().max() < 1e-12

    def test_rcbd_residual_df_and_mean(self):
        d = ag.residual_diagnostics(load_rcbd_example(), "yield", "treatment", block="block")
        assert abs(d.residuals["residual"].mean()) < 1e-12
        # sum of squares equals ANOVA error SS
        a = ag.anova_rcbd(load_rcbd_example(), "yield", "treatment", "block")
        sse = a.table.set_index("source").loc["Error", "ss"]
        assert (d.residuals["residual"] ** 2).sum() == pytest.approx(sse)

    def test_statistics_in_range_and_notes(self):
        d = ag.residual_diagnostics(load_crd_example(), "yield", "treatment")
        assert 0 <= d.shapiro_p <= 1
        assert 0 <= d.levene_p <= 1
        assert "does not prove" in d.notes
        assert "Shapiro" in str(d)

    def test_input_not_mutated(self):
        df = load_crd_example()
        before = df.copy()
        ag.residual_diagnostics(df, "yield", "treatment")
        pd.testing.assert_frame_equal(df, before)

    def test_errors(self):
        df = load_crd_example()
        with pytest.raises(AnalysisError):
            ag.residual_diagnostics(df, "nope", "treatment")
        with pytest.raises(AnalysisError):
            ag.residual_diagnostics(df[df["treatment"] == "N50"], "yield", "treatment")
        with pytest.raises(AnalysisError):
            ag.residual_diagnostics(load_rcbd_example().iloc[1:], "yield", "treatment", block="block")


class TestPlots:
    @pytest.mark.parametrize(
        "design",
        [
            lambda: ag.design.crd(n_treatments=3, replications=4, seed=1),
            lambda: ag.design.crd(n_treatments=3, replications=4, seed=1, row_length=4),
            lambda: ag.design.rcbd(n_treatments=4, replications=3, seed=1),
            lambda: ag.design.latin_square(n_treatments=4, seed=1),
            lambda: ag.design.factorial({"A": ["a", "b"], "B": ["x", "y"]}, replications=2, layout="rcbd"),
            lambda: ag.design.split_plot(["a", "b"], ["x", "y", "z"], replications=2),
        ],
    )
    def test_plot_layout_draws_every_plot(self, design, tmp_path):
        d = design()
        fig, ax = d.plot_layout()
        texts = {t.get_text().split("\n")[0] for t in ax.texts}
        assert set(d.table["plot_id"]) <= texts
        out = tmp_path / "layout.png"
        fig.savefig(out)
        assert out.stat().st_size > 0
        svg = tmp_path / "layout.svg"
        fig.savefig(svg)
        assert svg.stat().st_size > 0

    def test_plot_layout_does_not_change_table(self):
        d = ag.design.rcbd(n_treatments=3, replications=2, seed=1)
        before = d.table
        d.plot_layout()
        pd.testing.assert_frame_equal(before, d.table)

    def test_means_boxplot_interaction_residuals(self):
        df = load_rcbd_example()
        for err in ("se", "sd", "ci95"):
            fig, ax = ag.visualization.plot_means(df, "yield", "treatment", error=err)
            assert err.upper()[:2] in ax.get_ylabel().upper() or "CI" in ax.get_ylabel()
        letters = ag.tukey_hsd(df, "yield", "treatment", block="block").letters()
        ag.visualization.plot_means(df, "yield", "treatment", letters=letters)
        ag.visualization.plot_boxplot(df, "yield", "treatment")
        ag.visualization.plot_interaction(df, "yield", "treatment", "block")
        diag = ag.residual_diagnostics(df, "yield", "treatment", block="block")
        fig, axes = ag.visualization.plot_residuals(diag)
        assert len(axes) == 2

    def test_bad_inputs(self):
        df = load_rcbd_example()
        with pytest.raises(AnalysisError):
            ag.visualization.plot_means(df, "yield", "treatment", error="bad")
        with pytest.raises(AnalysisError):
            ag.visualization.plot_boxplot(df, "nope", "treatment")


class TestReports:
    @pytest.fixture
    def parts(self):
        df = load_rcbd_example()
        return {
            "design": ag.design.rcbd(n_treatments=4, replications=3, seed=1),
            "summary": ag.treatment_summary(df, "yield", "treatment"),
            "anova": ag.anova_rcbd(df, "yield", "treatment", "block"),
            "comparison": ag.tukey_hsd(df, "yield", "treatment", block="block"),
            "diagnostics": ag.residual_diagnostics(df, "yield", "treatment", block="block"),
        }

    def test_markdown(self, tmp_path, parts):
        p = ag.io.write_report(tmp_path / "r.md", title="Trial 1", **parts)
        text = p.read_text(encoding="utf-8")
        assert text.startswith("# Trial 1")
        for h in ("## Design", "## ANOVA", "Tukey HSD", "Compact letter display", "Residual diagnostics"):
            assert h in text
        assert "| source |" in text

    def test_html_escapes_title(self, tmp_path, parts):
        p = ag.io.write_report(tmp_path / "r.html", title="<b>x</b>", anova=parts["anova"])
        text = p.read_text(encoding="utf-8")
        assert "&lt;b&gt;x&lt;/b&gt;" in text
        assert "<b>x</b>" not in text
        assert "<table" in text

    def test_errors(self, tmp_path):
        with pytest.raises(ExportError):
            ag.io.write_report(tmp_path / "r.md")
        with pytest.raises(ExportError):
            ag.io.write_report(tmp_path / "r.pdf", anova=object())


class TestExtraDesignOptions:
    def test_crd_row_column(self):
        d = ag.design.crd(n_treatments=3, replications=4, seed=1, row_length=4)
        t = d.table
        assert t["row"].max() == 3
        assert t["column"].max() == 4
        assert not t.duplicated(["row", "column"]).any()
        assert "row" not in ag.design.crd(n_treatments=3, replications=4, seed=1).table.columns

    @pytest.mark.parametrize("bad", [0, -1, True, 2.5])
    def test_crd_bad_row_length(self, bad):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=2, replications=2, row_length=bad)

    @pytest.mark.parametrize("bad", [-3, True, "1", 1.5])
    def test_bad_start_plot(self, bad):
        with pytest.raises(DesignError):
            ag.design.rcbd(n_treatments=2, replications=2, start_plot=bad)

    def test_bad_prefix(self):
        with pytest.raises(DesignError):
            ag.design.crd(n_treatments=2, replications=2, plot_prefix=5)

    def test_table_mode_layout(self):
        d = ag.design.rcbd(n_treatments=3, replications=2, seed=1)
        out = d.text_layout(mode="table")
        assert "plot_id" in out
        assert "position_in_block" in out

    def test_bad_row_length_in_layout(self):
        d = ag.design.crd(n_treatments=2, replications=2, seed=1)
        with pytest.raises(ValueError, match="row_length"):
            d.text_layout(row_length=0)

    def test_rcbd_separator_character_used(self):
        d = ag.design.rcbd(n_treatments=3, replications=2, seed=1)
        assert any(set(line) == {"#"} for line in d.text_layout(separator="#").splitlines() if line)

    def test_incomplete_rcbd_anova_rejected(self):
        df = load_rcbd_example().iloc[1:]
        with pytest.raises(AnalysisError):
            ag.anova_rcbd(df, "yield", "treatment", "block")


class TestExcelExport:
    @pytest.fixture
    def parts(self):
        df = load_rcbd_example()
        return {
            "design": ag.design.rcbd(n_treatments=4, replications=3, seed=1),
            "summary": ag.treatment_summary(df, "yield", "treatment"),
            "anova": ag.anova_rcbd(df, "yield", "treatment", "block"),
            "comparison": ag.tukey_hsd(df, "yield", "treatment", block="block"),
            "diagnostics": ag.residual_diagnostics(df, "yield", "treatment", block="block"),
        }

    def test_all_sheets_and_values(self, tmp_path, parts):
        openpyxl = pytest.importorskip("openpyxl")
        p = ag.io.write_excel(tmp_path / "out.xlsx", **parts)
        wb = openpyxl.load_workbook(p)
        assert wb.sheetnames == [
            "Design", "Layout", "Summary", "ANOVA", "ANOVA notes",
            "Comparison", "Letters", "Comparison notes", "Residuals", "Diagnostics",
        ]
        back = pd.read_excel(p, sheet_name="ANOVA")
        pd.testing.assert_frame_equal(back, parts["anova"].table, check_dtype=False)
        assert len(pd.read_excel(p, sheet_name="Design")) == 12

    def test_report_xlsx_delegates(self, tmp_path, parts):
        pytest.importorskip("openpyxl")
        p = ag.io.write_report(tmp_path / "r.xlsx", anova=parts["anova"])
        assert pd.ExcelFile(p).sheet_names == ["ANOVA", "ANOVA notes"]

    def test_nothing_supplied(self, tmp_path):
        with pytest.raises(ExportError):
            ag.io.write_excel(tmp_path / "x.xlsx")

    def test_unwritable_path(self, tmp_path, parts):
        pytest.importorskip("openpyxl")
        with pytest.raises(ExportError):
            ag.io.write_excel(tmp_path / "missing_dir" / "x.xlsx", anova=parts["anova"])

"""Tests for the one-command analyze() workflow."""
import pandas as pd
import pytest

import agristatspy as ag
from agristatspy.datasets import (
    load_crd_example,
    load_factorial_example,
    load_latin_square_example,
    load_rcbd_example,
    load_split_plot_example,
)
from agristatspy.exceptions import AnalysisError, ExportError

NON_NORMAL = pd.DataFrame(
    {
        "t": list("A" * 6 + "B" * 6 + "C" * 6),
        "y": [5.0, 5.2, 4.9, 5.1, 5.0, 12.0, 6.0, 6.1, 5.9, 6.2, 6.0, 6.1,
              7.0, 7.1, 6.9, 7.2, 7.0, 7.1],
    }
)


@pytest.fixture
def crd_report():
    return ag.analyze(load_crd_example(), "yield", "crd", treatment="treatment")


class TestNormalPath:
    def test_runs_everything_in_one_call(self, crd_report):
        r = crd_report
        assert r.status == "complete"
        assert r.normal is True
        assert r.anova is not None
        assert r.comparison is not None
        assert list(r.means["group"]) == ["c", "b", "a"]  # Control, N50, N100

    def test_stats_match_separate_calls(self, crd_report):
        data = load_crd_example()
        a = ag.anova_crd(data, "yield", "treatment")
        assert crd_report.stats["grand_mean"] == pytest.approx(a.grand_mean)
        assert crd_report.stats["cv_pct"] == pytest.approx(a.cv_pct)
        mse = a.table.set_index("source").loc["Error", "ms"]
        assert crd_report.stats["sem_pooled"] == pytest.approx((mse / 4) ** 0.5)

    def test_sem_and_cd_hand_check(self, crd_report):
        m = crd_report.means.set_index("treatment")
        ctrl = load_crd_example().query("treatment == 'Control'")["yield"]
        assert m.loc["Control", "mean"] == pytest.approx(ctrl.mean())
        assert m.loc["Control", "sem"] == pytest.approx(ctrl.std(ddof=1) / 2)
        from scipy import stats
        s = crd_report.stats
        expected_cd = stats.t.ppf(0.975, s["df_error"]) * (2 * s["mse"] / 4) ** 0.5
        assert s["cd"] == pytest.approx(expected_cd)

    def test_rcbd_uses_block_error_term(self):
        r = ag.analyze(load_rcbd_example(), "yield", "rcbd", treatment="treatment", block="block",
                       on_nonnormal="proceed")
        assert r.comparison.df_error == 6
        assert r.status == "complete"

    def test_latin_square(self):
        r = ag.analyze(load_latin_square_example(), "yield", "latin_square",
                       treatment="treatment", row="row", column="column")
        assert r.status == "complete"
        assert r.comparison.df_error == 6

    def test_factorial_cells(self):
        r = ag.analyze(load_factorial_example(), "yield", "factorial",
                       factors=["variety", "nitrogen"], block="block", on_nonnormal="proceed")
        assert len(r.means) == 6
        assert "group" in r.means
        assert "ANOVA p (variety:nitrogen)" in set(r.thesis_table()["Treatment"])

    def test_split_plot_no_posthoc(self):
        r = ag.analyze(load_split_plot_example(), "yield", "split_plot", main="irrigation",
                       sub="variety", block="block", on_nonnormal="proceed")
        assert r.comparison is None
        assert "group" not in r.means
        assert any("split-plot" in m for m in r.messages)

    def test_posthoc_skipped_when_not_significant(self):
        df = pd.DataFrame({"t": list("ABC") * 4,
                           "y": [5.0, 5.1, 4.9, 5.1, 4.9, 5.0, 5.0, 5.05, 4.95, 4.95, 5.05, 5.0]})
        r = ag.analyze(df, "y", "crd", treatment="t", on_nonnormal="proceed")
        assert r.comparison is None
        assert any("not run" in m for m in r.messages)
        r2 = ag.analyze(df, "y", "crd", treatment="t", on_nonnormal="proceed", posthoc_when="always")
        assert r2.comparison is not None

    @pytest.mark.parametrize("method", ["lsd", "bonferroni", None])
    def test_posthoc_choice(self, method):
        r = ag.analyze(load_crd_example(), "yield", "crd", treatment="treatment", posthoc=method)
        assert (r.comparison is None) == (method is None)
        if method:
            assert method[:3] in r.comparison.method.lower()


class TestNonNormalDecision:
    def test_non_interactive_pauses_then_proceeds(self):
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t", ask=None, on_nonnormal="ask")
        # pytest stdin is not a tty -> paused
        assert r.status == "paused"
        assert r.normal is False
        assert r.anova is None
        assert len(r.flagged) == 1
        assert r.flagged.iloc[0]["y"] == 12.0
        done = r.proceed()
        assert done.status == "complete"
        assert done.anova is not None
        assert r.status == "paused"  # original untouched

    def test_ask_proceed(self):
        asked = []
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t",
                       ask=lambda q: asked.append(q) or "p")
        assert r.status == "complete"
        assert "NOT normal" in asked[0]
        assert any("proceeding" in m for m in r.messages)

    def test_ask_review(self):
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t", ask=lambda q: "r")
        assert r.status == "review"
        assert r.anova is None
        assert len(r.flagged) == 1
        assert "transformation" in " ".join(r.messages)

    def test_ask_reprompts_on_invalid_answer(self):
        answers = iter(["maybe", "", "p"])
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t", ask=lambda q: next(answers))
        assert r.status == "complete"

    def test_flags_policy_values(self):
        assert ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="proceed").status == "complete"
        assert ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="review").status == "review"
        assert ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="stop").status == "stopped"

    def test_data_never_modified(self):
        before = NON_NORMAL.copy()
        ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="proceed")
        pd.testing.assert_frame_equal(NON_NORMAL, before)

    def test_tables_unavailable_until_proceed(self):
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="stop")
        with pytest.raises(AnalysisError):
            r.thesis_table()
        with pytest.raises(ExportError):
            r.export("x.csv")

    def test_str_mentions_status_and_flags(self):
        r = ag.analyze(NON_NORMAL, "y", "crd", treatment="t", on_nonnormal="review")
        text = str(r)
        assert "NOT normal" in text
        assert "kept, not removed" in text


class TestTablesAndExport:
    def test_thesis_table(self, crd_report):
        t = crd_report.thesis_table(decimals=2)
        assert list(t.columns) == ["Treatment", "Mean ± SEM", "Group"]
        assert t.iloc[0]["Mean ± SEM"] == "2.48 ± 0.09"
        labels = list(t["Treatment"])
        assert "SEm ±" in labels
        assert "CD (5%)" in labels
        assert "CV (%)" in labels
        assert len(crd_report.thesis_table(footer=False)) == 3

    def test_anova_table_stars(self, crd_report):
        t = crd_report.anova_table()
        assert list(t.columns) == ["Source", "df", "SS", "MS", "F", "p", "Sig."]
        assert t.loc[t["Source"] == "Treatment", "Sig."].iloc[0] == "***"
        assert t.loc[t["Source"] == "Treatment", "p"].iloc[0] == "<0.001"

    @pytest.mark.parametrize("ext", ["csv", "md", "html", "tex"])
    def test_export_text_formats(self, tmp_path, crd_report, ext):
        p = crd_report.export(tmp_path / f"out.{ext}")
        text = p.read_text(encoding="utf-8")
        assert "N100" in text
        if ext == "tex":
            assert "\\begin{tabular}" in text
            assert "$\\pm$" in text
            assert "CD (5\\%)" in text
        if ext == "html":
            assert "<table" in text

    def test_export_xlsx(self, tmp_path, crd_report):
        pytest.importorskip("openpyxl")
        p = crd_report.export(tmp_path / "out.xlsx")
        sheets = pd.ExcelFile(p).sheet_names
        assert {"Treatment means", "ANOVA", "Pairwise", "Summary", "Report"} <= set(sheets)

    def test_bad_format(self, tmp_path, crd_report):
        with pytest.raises(ExportError):
            crd_report.export(tmp_path / "out.pdf")


class TestValidation:
    def test_unknown_design(self):
        with pytest.raises(AnalysisError):
            ag.analyze(load_crd_example(), "yield", "mystery", treatment="treatment")

    def test_missing_args(self):
        with pytest.raises(AnalysisError):
            ag.analyze(load_rcbd_example(), "yield", "rcbd", treatment="treatment")

    def test_missing_column(self):
        with pytest.raises(AnalysisError):
            ag.analyze(load_crd_example(), "nope", "crd", treatment="treatment")

    @pytest.mark.parametrize(
        "kw",
        [{"alpha": 0}, {"alpha": 1.5}, {"normality_alpha": True}, {"posthoc": "duncan"},
         {"posthoc_when": "sometimes"}, {"on_nonnormal": "maybe"}],
    )
    def test_bad_options(self, kw):
        with pytest.raises(AnalysisError):
            ag.analyze(load_crd_example(), "yield", "crd", treatment="treatment", **kw)

    def test_incomplete_rcbd_rejected(self):
        with pytest.raises(AnalysisError):
            ag.analyze(load_rcbd_example().iloc[1:], "yield", "rcbd", treatment="treatment", block="block")

import pytest

import agristatspy as ag
from agristatspy.datasets import load_rcbd_example
from agristatspy.exceptions import AnalysisError


def test_dmrt_basic_and_between_lsd_and_tukey():
    df = load_rcbd_example()
    kw = {"response": "yield", "treatment": "treatment", "block": "block"}
    n = {f.__name__: f(df, **kw).table["reject"].sum() for f in (ag.lsd_test, ag.dmrt, ag.tukey_hsd)}
    assert n["lsd_test"] >= n["dmrt"] >= n["tukey_hsd"]
    res = ag.dmrt(df, **kw)
    assert res.method == "Duncan DMRT"
    assert set(res.letters().columns) == {"treatment", "mean", "group"}


def test_analyze_all_methods():
    df = load_rcbd_example()
    rep = ag.analyze(df, "yield", "rcbd", treatment="treatment", block="block",
                     posthoc="all", posthoc_when="always", on_nonnormal="proceed")
    assert set(rep.comparisons) == {"tukey", "lsd", "dmrt", "bonferroni"}
    t = rep.thesis_table()
    assert sum(c.startswith("Group") for c in t.columns) == 4
    assert "Duncan DMRT" in str(rep)


def test_analyze_list_and_ask():
    df = load_rcbd_example()
    kw = {"treatment": "treatment", "block": "block", "posthoc_when": "always", "on_nonnormal": "proceed"}
    rep = ag.analyze(df, "yield", "rcbd", posthoc=["lsd", "dmrt"], **kw)
    assert list(rep.comparisons) == ["lsd", "dmrt"]
    rep = ag.analyze(df, "yield", "rcbd", posthoc="ask", ask=lambda q: "3", **kw)
    assert list(rep.comparisons) == ["dmrt"]
    rep = ag.analyze(df, "yield", "rcbd", posthoc="ask", ask=lambda q: "lsd,dmrt", **kw)
    assert list(rep.comparisons) == ["lsd", "dmrt"]


def test_bad_posthoc():
    df = load_rcbd_example()
    with pytest.raises(AnalysisError):
        ag.analyze(df, "yield", "rcbd", treatment="treatment", block="block", posthoc="scheffe")

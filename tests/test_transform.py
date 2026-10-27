"""Tests for ag.transform, checked against hand-computed values."""
import math

import numpy as np
import pandas as pd
import pytest

import agristatspy as ag
from agristatspy.exceptions import AnalysisError


def df(values):
    return pd.DataFrame({"y": values})


def test_log_hand_values_and_original_kept():
    d = df([1.0, math.e, math.e**2])
    r = ag.transform(d, "y", "log")
    assert r.column == "y_log"
    np.testing.assert_allclose(r.data["y_log"], [0.0, 1.0, 2.0])
    assert list(d.columns) == ["y"]  # input untouched
    assert list(r.data["y"]) == list(d["y"])


def test_log10_with_constant():
    r = ag.transform(df([0, 9, 99]), "y", "log10", constant=1)
    np.testing.assert_allclose(r.data["y_log10"], [0.0, 1.0, 2.0])


def test_sqrt_counts():
    r = ag.transform(df([0, 3.5]), "y", "sqrt", constant=0.5)
    np.testing.assert_allclose(r.data["y_sqrt"], [math.sqrt(0.5), 2.0])


def test_arcsine_percent_and_proportion():
    r = ag.transform(df([0, 25, 50, 100]), "y", "arcsine")
    np.testing.assert_allclose(r.data["y_arcsine"], [0, 30, 45, 90], atol=1e-9)
    p = ag.transform(df([0.25]), "y", "arcsine", scale="proportion")
    assert p.data["y_arcsine"][0] == pytest.approx(30.0)


def test_reciprocal_and_square():
    np.testing.assert_allclose(
        ag.transform(df([2, 4]), "y", "reciprocal").data["y_reciprocal"], [0.5, 0.25]
    )
    np.testing.assert_allclose(ag.transform(df([2, 3]), "y", "square").data["y_square"], [4, 9])


def test_boxcox_fixed_lambda_matches_formula():
    r = ag.transform(df([1.0, 2.0, 4.0]), "y", "boxcox", lmbda=0.5)
    np.testing.assert_allclose(r.data["y_boxcox"], [(x**0.5 - 1) / 0.5 for x in (1, 2, 4)])
    assert r.params["lambda"] == 0.5


def test_boxcox_estimated_lambda_matches_scipy():
    from scipy import stats

    x = np.array([1.2, 2.5, 3.1, 7.9, 15.0, 22.0, 40.0])
    r = ag.transform(df(x), "y", "boxcox")
    expected, lam = stats.boxcox(x)
    assert r.params["lambda"] == pytest.approx(lam)
    np.testing.assert_allclose(r.data["y_boxcox"], expected)


@pytest.mark.parametrize(
    ("method", "kw", "values"),
    [
        ("log", {}, [0.5, 3.0, 12.0]),
        ("log", {"constant": 1.0}, [0.0, 3.0]),
        ("log10", {}, [0.5, 3.0]),
        ("sqrt", {"constant": 0.5}, [0.0, 3.0]),
        ("arcsine", {}, [5.0, 40.0, 95.0]),
        ("arcsine", {"scale": "proportion"}, [0.05, 0.4]),
        ("reciprocal", {}, [0.5, 3.0]),
        ("square", {}, [0.5, 3.0]),
        ("boxcox", {"lmbda": 0.0}, [0.5, 3.0]),
        ("boxcox", {"lmbda": 0.5}, [0.5, 3.0]),
        ("boxcox", {}, [0.5, 3.0, 8.0, 20.0]),
    ],
)
def test_back_transform_round_trip(method, kw, values):
    r = ag.transform(df(values), "y", method, **kw)
    np.testing.assert_allclose(r.back_transform(r.data[r.column]), values, atol=1e-9)
    assert isinstance(r.back_transform(float(r.data[r.column].iloc[0])), float)


def test_missing_values_preserved():
    r = ag.transform(df([1.0, np.nan, 4.0]), "y", "sqrt")
    assert np.isnan(r.data["y_sqrt"][1])
    assert len(r.data) == 3


@pytest.mark.parametrize(
    ("method", "kw", "values"),
    [
        ("log", {}, [0.0, 1.0]),
        ("log", {}, [-1.0, 1.0]),
        ("sqrt", {}, [-0.1, 1.0]),
        ("reciprocal", {}, [0.0, 1.0]),
        ("arcsine", {}, [101.0]),
        ("arcsine", {"scale": "proportion"}, [20.0]),
        ("boxcox", {}, [0.0, 1.0]),
    ],
)
def test_out_of_domain_raises(method, kw, values):
    with pytest.raises(AnalysisError):
        ag.transform(df(values), "y", method, **kw)


def test_bad_inputs_raise():
    d = pd.DataFrame({"y": [1.0, 2.0], "t": ["a", "b"], "y_log": [0, 0]})
    with pytest.raises(AnalysisError):
        ag.transform(d, "y", "nope")
    with pytest.raises(AnalysisError):
        ag.transform(d, "missing", "log")
    with pytest.raises(AnalysisError):
        ag.transform(d, "t", "log")
    with pytest.raises(AnalysisError):
        ag.transform(d, "y", "log")  # y_log exists
    with pytest.raises(AnalysisError):
        ag.transform(d, "y", "arcsine", scale="bad")
    assert ag.transform(d, "y", "log", name="ly").column == "ly"


def test_feeds_into_anova():
    d = pd.DataFrame(
        {"trt": ["A"] * 4 + ["B"] * 4, "y": [1, 2, 1.5, 2.5, 10, 20, 15, 25.0]}
    )
    r = ag.transform(d, "y", "log")
    res = ag.anova_crd(r.data, response=r.column, treatment="trt")
    assert res is not None
    assert "log" in str(r)

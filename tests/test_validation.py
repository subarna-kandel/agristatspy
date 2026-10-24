"""
Tests for input validation utilities (shared by CRD and RCBD).
"""
import pytest

from agristatspy.design.validation import (
    parse_treatments,
    validate_replications,
    validate_seed,
)
from agristatspy.exceptions import DesignError


class TestParseTreatments:
    def test_explicit_labels(self):
        result = parse_treatments(["A", "B", "C"], None)
        assert result == ["A", "B", "C"]

    def test_auto_labels(self):
        result = parse_treatments(None, 4)
        assert result == ["T1", "T2", "T3", "T4"]

    def test_both_consistent(self):
        result = parse_treatments(["A", "B", "C"], 3)
        assert result == ["A", "B", "C"]

    def test_both_inconsistent_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(["A", "B", "C"], 2)

    def test_neither_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(None, None)

    def test_empty_list_raises(self):
        with pytest.raises(DesignError):
            parse_treatments([], None)

    def test_empty_label_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(["A", "  ", "B"], None)

    def test_duplicate_labels_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(["A", "B", "A"], None)

    def test_zero_n_treatments_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(None, 0)

    def test_bool_n_treatments_raises(self):
        with pytest.raises(DesignError):
            parse_treatments(None, True)


class TestValidateReplications:
    def test_valid_integer(self):
        assert validate_replications(3) == 3

    def test_zero_raises(self):
        with pytest.raises(DesignError):
            validate_replications(0)

    def test_negative_raises(self):
        with pytest.raises(DesignError):
            validate_replications(-5)

    def test_bool_raises(self):
        with pytest.raises(DesignError):
            validate_replications(True)

    def test_float_raises(self):
        with pytest.raises(DesignError):
            validate_replications(2.0)


class TestValidateSeed:
    def test_none_passes(self):
        assert validate_seed(None) is None

    def test_integer_passes(self):
        assert validate_seed(42) == 42

    def test_zero_passes(self):
        assert validate_seed(0) == 0

    def test_bool_raises(self):
        with pytest.raises(DesignError):
            validate_seed(True)

    def test_float_raises(self):
        with pytest.raises(DesignError):
            validate_seed(3.14)

    def test_string_raises(self):
        with pytest.raises(DesignError):
            validate_seed("42")

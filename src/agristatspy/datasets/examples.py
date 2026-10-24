"""
Built-in example datasets for testing and tutorials.

These are small hand-crafted datasets designed to produce clean, verifiable
ANOVA results.  They are not intended to represent real field-trial data.
"""
from __future__ import annotations

import pandas as pd


def load_crd_example() -> pd.DataFrame:
    """Return a small CRD dataset with 3 treatments × 4 replications.

    Response variable: ``yield`` (simulated grain yield, t/ha).
    The data are chosen so that the treatment effect is detectable.

    Treatment means (approximate):
    - Control ≈ 2.5
    - N50     ≈ 3.5
    - N100    ≈ 4.5

    Returns
    -------
    pd.DataFrame
        Columns: ``plot_id``, ``treatment``, ``replication_number``, ``yield``.
    """
    records = [
        # Control
        {"plot_id": "P01", "treatment": "Control", "replication_number": 1, "yield": 2.3},
        {"plot_id": "P02", "treatment": "Control", "replication_number": 2, "yield": 2.7},
        {"plot_id": "P03", "treatment": "Control", "replication_number": 3, "yield": 2.5},
        {"plot_id": "P04", "treatment": "Control", "replication_number": 4, "yield": 2.4},
        # N50
        {"plot_id": "P05", "treatment": "N50", "replication_number": 1, "yield": 3.4},
        {"plot_id": "P06", "treatment": "N50", "replication_number": 2, "yield": 3.6},
        {"plot_id": "P07", "treatment": "N50", "replication_number": 3, "yield": 3.5},
        {"plot_id": "P08", "treatment": "N50", "replication_number": 4, "yield": 3.3},
        # N100
        {"plot_id": "P09", "treatment": "N100", "replication_number": 1, "yield": 4.4},
        {"plot_id": "P10", "treatment": "N100", "replication_number": 2, "yield": 4.6},
        {"plot_id": "P11", "treatment": "N100", "replication_number": 3, "yield": 4.5},
        {"plot_id": "P12", "treatment": "N100", "replication_number": 4, "yield": 4.5},
    ]
    return pd.DataFrame(records)


def load_rcbd_example() -> pd.DataFrame:
    """Return a small RCBD dataset with 4 treatments × 3 blocks.

    Response variable: ``yield`` (simulated grain yield, t/ha).

    Returns
    -------
    pd.DataFrame
        Columns: ``plot_id``, ``block``, ``treatment``, ``yield``.
    """
    records = [
        # Block 1
        {"plot_id": "P01", "block": 1, "treatment": "Control", "yield": 2.1},
        {"plot_id": "P02", "block": 1, "treatment": "N50",     "yield": 3.2},
        {"plot_id": "P03", "block": 1, "treatment": "N100",    "yield": 4.1},
        {"plot_id": "P04", "block": 1, "treatment": "N150",    "yield": 4.8},
        # Block 2
        {"plot_id": "P05", "block": 2, "treatment": "Control", "yield": 2.3},
        {"plot_id": "P06", "block": 2, "treatment": "N50",     "yield": 3.5},
        {"plot_id": "P07", "block": 2, "treatment": "N100",    "yield": 4.3},
        {"plot_id": "P08", "block": 2, "treatment": "N150",    "yield": 5.0},
        # Block 3
        {"plot_id": "P09", "block": 3, "treatment": "Control", "yield": 2.0},
        {"plot_id": "P10", "block": 3, "treatment": "N50",     "yield": 3.1},
        {"plot_id": "P11", "block": 3, "treatment": "N100",    "yield": 4.0},
        {"plot_id": "P12", "block": 3, "treatment": "N150",    "yield": 4.7},
    ]
    return pd.DataFrame(records)

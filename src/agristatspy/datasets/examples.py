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


# Fixed, hand-picked "noise" values so the example datasets are deterministic.
_NOISE = [0.12, -0.08, 0.05, -0.15, 0.10, -0.02, 0.07, -0.11, 0.03, 0.09,
          -0.06, 0.14, -0.13, 0.01, 0.08, -0.04, 0.11, -0.09, 0.02, -0.07,
          0.06, -0.01, 0.13, -0.12, 0.04, -0.05, 0.09, -0.10, 0.00, 0.15,
          -0.03, 0.07, -0.14, 0.05, 0.12, -0.08, 0.01, 0.10, -0.06, 0.03,
          -0.11, 0.08, -0.02, 0.14, -0.09, 0.04, -0.13, 0.06, 0.02, -0.07]


def load_latin_square_example() -> pd.DataFrame:
    """Return a 4 × 4 Latin square dataset (fertilizer on a sloping field).

    Rows model the slope direction and columns the cross-slope direction.

    Returns
    -------
    pd.DataFrame
        Columns: ``plot_id``, ``row``, ``column``, ``treatment``, ``yield``.
    """
    trts = ["Control", "N50", "N100", "N150"]
    effect = {"Control": 2.0, "N50": 3.0, "N100": 3.8, "N150": 4.2}
    records = []
    for r in range(4):
        for c in range(4):
            trt = trts[(r + c) % 4]
            n = _NOISE[r * 4 + c]
            records.append(
                {
                    "plot_id": f"P{r * 4 + c + 1:02d}",
                    "row": r + 1,
                    "column": c + 1,
                    "treatment": trt,
                    "yield": round(effect[trt] + 0.15 * r - 0.1 * c + n, 3),
                }
            )
    return pd.DataFrame(records)


def load_factorial_example() -> pd.DataFrame:
    """Return a 2 × 3 factorial (variety × nitrogen) RCBD dataset, 3 blocks.

    Returns
    -------
    pd.DataFrame
        Columns: ``plot_id``, ``block``, ``variety``, ``nitrogen``,
        ``treatment``, ``yield``.
    """
    records = []
    k = 0
    for b in (1, 2, 3):
        for v, v_eff in (("V1", 0.0), ("V2", 0.6)):
            for n, n_eff in (("N0", 0.0), ("N50", 1.0), ("N100", 1.5)):
                inter = 0.3 if (v == "V2" and n == "N100") else 0.0
                records.append(
                    {
                        "plot_id": f"P{k + 1:02d}",
                        "block": b,
                        "variety": v,
                        "nitrogen": n,
                        "treatment": f"{v} | {n}",
                        "yield": round(3.0 + v_eff + n_eff + inter + 0.2 * b + _NOISE[k], 3),
                    }
                )
                k += 1
    return pd.DataFrame(records)


def load_split_plot_example() -> pd.DataFrame:
    """Return a split-plot dataset: irrigation (main) × variety (sub), 3 blocks.

    Returns
    -------
    pd.DataFrame
        Columns: ``plot_id``, ``block``, ``main_plot``, ``irrigation``,
        ``variety``, ``yield``.
    """
    records = []
    k = 0
    for b in (1, 2, 3):
        for m, (irr, i_eff) in enumerate((("Rainfed", 0.0), ("Drip", 1.2)), start=1):
            whole_plot_noise = _NOISE[(b * 7 + m * 3) % len(_NOISE)]
            for var, v_eff in (("V1", 0.0), ("V2", 0.5), ("V3", 0.9)):
                records.append(
                    {
                        "plot_id": f"P{k + 1:02d}",
                        "block": b,
                        "main_plot": f"B{b}-M{m}",
                        "irrigation": irr,
                        "variety": var,
                        "yield": round(4.0 + i_eff + v_eff + 0.25 * b + whole_plot_noise + _NOISE[k], 3),
                    }
                )
                k += 1
    return pd.DataFrame(records)

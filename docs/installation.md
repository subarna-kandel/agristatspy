# Installation Guide

## Requirements

- Python 3.10, 3.11, or 3.12
- Standard scientific stack:
  - `pandas >= 2.0`
  - `numpy >= 1.24`
  - `scipy >= 1.10`
  - `statsmodels >= 0.14`
  - `matplotlib >= 3.7`

---

## Installing AgriStatsPy

### Standard Installation

Install via `pip`:

```bash
pip install agristatspy
```

### With Excel Support

To export experimental plans directly to `.xlsx` files:

```bash
pip install agristatspy[excel]
```

### Development Installation

Clone the repository and install with development dependencies:

```bash
git clone https://github.com/agristatspy/agristatspy.git
cd agristatspy
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Run tests to verify the setup:

```bash
pytest
```

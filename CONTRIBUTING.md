# Contributing to AgriStatsPy

Thanks for your interest! Please open an issue before starting work on a new
statistical feature so the method can be agreed first.

## Setup

```bash
git clone https://github.com/agristatspy/agristatspy.git
cd agristatspy
pip install -e ".[dev]"
pytest
ruff check src tests examples
```

## Requirements for new features

1. A written method: model, assumptions, limitations, and error behaviour.
2. Tests against hand-checkable values or an independent reference
   (`scipy`, `statsmodels`, R, or published examples).
3. Input validation with clear `AgriStatsError` subclasses; never silently
   drop data or change requested counts.
4. Documentation (`docs/`), a `CHANGELOG.md` entry, and an example if useful.
5. Statistical output must state the model, test, alpha, and missing-data
   handling, and must not overstate what significance or diagnostics prove.

## Pull requests

Keep changes focused, make sure `pytest` and `ruff` pass, and describe the
statistical reasoning for any numerical change.

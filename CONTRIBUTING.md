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

## Guidelines for AI agents and AI-assisted contributions

AI-assisted contributions are welcome, but they are held to the same standard
as any other and must follow these strict rules. Agents (Claude Code, Copilot,
Cursor, etc.) must read this section before making changes.

### Mandatory rules

1. **Disclose AI involvement.** Every pull request and every related issue
   must state that AI was used, which tool/model, and exactly what it did
   (e.g. "implemented `foo()`, wrote tests, updated docs"). Do not hide or
   understate AI work.
2. **Link the issue.** Every PR must reference an issue (`Closes #123`). Do not
   start work on a new statistical method before the method is agreed in an
   issue.
3. **Run everything before you push or open a PR.** Do not commit, push, or
   open a PR until the code has actually been run and verified locally:
   - `pip install -e ".[dev]"`
   - `pytest` passes in full
   - `ruff check src tests examples` is clean
   - any new or changed example/script executes end to end
4. **Report results honestly.** Paste the real test/lint output (or a faithful
   summary) in the PR. If something fails or was not run, say so. Never claim
   tests pass without running them.
5. **Never fabricate numbers.** Reference values in tests must come from an
   independent source (`scipy`, `statsmodels`, R, published tables) or a
   hand calculation, and the source must be cited. Do not generate expected
   values from the code under test.
6. **No silent behaviour changes.** Do not alter numerical output, defaults,
   or public API without calling it out in the PR and `CHANGELOG.md`.
7. **Keep scope tight.** One focused change per PR. No unrelated refactors,
   reformatting, dependency additions, or deleted tests.
8. **A human is accountable.** A human contributor must review, understand,
   and take responsibility for every line submitted.

### PR checklist for AI-assisted work

```text
- [ ] AI use disclosed (tool/model and what it did)
- [ ] Linked issue (Closes #...)
- [ ] pytest passes locally (output included)
- [ ] ruff check passes locally
- [ ] Examples/scripts run end to end
- [ ] Tests use independent reference values (sources cited)
- [ ] Docs and CHANGELOG.md updated
- [ ] Reviewed line by line by a human
```

PRs that omit AI disclosure, have no linked issue, or were pushed without
being run may be closed without review.

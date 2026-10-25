# Changelog

All notable changes to this project are documented here. The project follows
[Semantic Versioning](https://semver.org/) and the
[Keep a Changelog](https://keepachangelog.com/) format.

## [Unreleased]

### Added
- `ag.analyze` is now the single entry point: `design="alpha_lattice"` (+ `replicate=`), `environment=` and a list of `response`s dispatch to the matching analysis; all honour `posthoc=` and `.export()`; `genetics=True` adds genetic parameters to crd/rcbd analyses.
- `ag.design.alpha_lattice` (cyclic alpha(0,1) generating array, randomized) with text and graphical layouts, and `ag.anova_alpha_lattice` (intra-block ANOVA, adjusted means, SED matrix, relative efficiency, `.compare()`).
- `ag.multi_environment` — combined ANOVA over locations/years/location × year (random or fixed environments), Bartlett homogeneity test, per-environment ANOVAs, `.compare()`, export.
- `ag.analyze_traits`, `ag.genetic_parameters` (GCV, PCV, ECV, h², GA, GAM, classes; single or multiple environments), `ag.correlation`, `ag.pca` (+ `plot_correlation`, `plot_pca`), `ag.two_sample_test`, `ag.compare_means`.
- `ag.dmrt` — Duncan's Multiple Range Test (CRD, RCBD, Latin square error terms).
- `ag.analyze(posthoc=...)` now accepts `"dmrt"`, a list of tests, `"all"` or
  `"ask"`; every test run is kept in `report.comparisons` and shown side by side
  in the thesis table and exports.

## [0.2.0] - 2026-10-24

### Added
- Designs: `latin_square`, `factorial` (CRD or RCBD layout), `split_plot`.
- `crd(row_length=...)` adds row/column positions; `text_layout(mode="table")`.
- ANOVA: `anova_latin_square`, `anova_factorial`, `anova_split_plot`, and the
  `anova(...)` dispatcher.
- Mean comparisons: `lsd_test`, `bonferroni`, compact letter displays
  (`ComparisonResult.letters()`), and `block=` support (RCBD or Latin square
  error terms) for all procedures.
- `residual_diagnostics` (Shapiro–Wilk, Brown–Forsythe).
- Graphics: `design.plot_layout()`, `plot_means`, `plot_boxplot`,
  `plot_interaction`, `plot_residuals` (PNG/SVG).
- Markdown/HTML reports via `ag.io.write_report`, and multi-sheet Excel
  workbooks via `ag.io.write_excel` (also `write_report(..., "x.xlsx")`).
- New example datasets and example scripts 05–08.
- `CHANGELOG.md`, `CONTRIBUTING.md`, `CITATION.cff`, MIT `LICENSE`, issue templates.

- `ag.analyze` one-command workflow (Shapiro–Wilk → decision → ANOVA → post-hoc letters → mean ± SEM, SEm±, CD, CV%) with `AnalysisReport.thesis_table()`, `anova_table()`, `proceed()` and `export()` (xlsx, csv, md, html, tex).
- `residual_diagnostics` accepts a list of blocking columns.

### Changed
- `tukey_hsd` is now computed directly from the studentized range distribution
  (verified against `statsmodels`) and validates `alpha`.
- `anova_rcbd` rejects layouts with missing treatment–block combinations.
- `crd`/`rcbd` validate `plot_prefix` and `start_plot`; text layouts validate
  `row_length`, honour `separator` for RCBD block dividers, and reject unknown
  RCBD modes.

## [0.1.0]

### Added
- CRD and RCBD generation, text layouts, CSV/Excel export, descriptive
  statistics, one-way and RCBD ANOVA, Tukey HSD, documentation and tests.

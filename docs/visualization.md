# Field Visualization Guide

AgriStatsPy provides human-readable visual layouts designed for field staff, agronomists, and researchers.

---

## 1. Text Layouts (Phase 1 MVP)

Field trials require unambiguous plot identification and clear borders between blocks. AgriStatsPy's `text_layout()` renders aligned ASCII/Unicode diagrams.

### RCBD Block Layouts with Dividers

When inspecting or printing RCBD plans, block separation is paramount:

```python
import agristatspy as ag

design = ag.design.rcbd(
    treatments=["Control", "N50", "N100", "N150"],
    replications=3,
    seed=42,
)

print(design.text_layout(block_separator=True, separator="="))
```

Output:
```text
RCBD — 4 treatments × 3 blocks = 12 plots
Seed: 42

BLOCK 1
+------------+------------+------------+------------+
| Plot P01   | Plot P02   | Plot P03   | Plot P04   |
| N50        | Control    | N150       | N100       |
+------------+------------+------------+------------+

=====================================================

BLOCK 2
+------------+------------+------------+------------+
| Plot P05   | Plot P06   | Plot P07   | Plot P08   |
| N100       | N50        | Control    | N150       |
+------------+------------+------------+------------+

=====================================================

BLOCK 3
+------------+------------+------------+------------+
| Plot P09   | Plot P10   | Plot P11   | Plot P12   |
| Control    | N150       | N100       | N50        |
+------------+------------+------------+------------+
```

### Customizing Cell Widths & Formatting

You can adjust `cell_width` or change the separator character:

```python
print(design.text_layout(cell_width=16, separator="-"))
```

---

## 2. Graphical Layouts

`design.plot_layout()` draws one cell per plot, coloured by treatment, with the plot ID and treatment inside each cell. RCBD, factorial-RCBD and split-plot designs get one outlined band per block; Latin squares are drawn as a row × column grid.

```python
fig, ax = design.plot_layout()
fig.savefig("field_layout.png", dpi=300, bbox_inches="tight")
fig.savefig("field_layout.svg")
```

Rendering is a view of the assignment table and never changes the randomization. A CRD grid is a numbering convention (use `row_length=` in `crd()` to store row/column positions); it does not model field gradients.

## 3. Result Plots

All functions return `(fig, ax)` and never call `plt.show()`.

| Function | Shows |
|---|---|
| `ag.visualization.plot_means(df, y, trt, error="se"/"sd"/"ci95", letters=None)` | Treatment means with the chosen error bars (named on the axis); optional letter groups |
| `ag.visualization.plot_boxplot(df, y, trt)` | Boxplot with the individual observations |
| `ag.visualization.plot_interaction(df, y, x, trace)` | Interaction (profile) plot |
| `ag.visualization.plot_residuals(diagnostics)` | Residuals vs fitted and normal Q–Q plot |

Interaction and residual plots are visual aids; interpret them with the cautions in `statistical-methods.md`.

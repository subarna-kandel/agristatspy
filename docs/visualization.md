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

## 2. Graphical Layouts (Phase 3)

Graphical plot maps using `matplotlib` (such as heatmaps, colored plot cards, and spatial field maps) are scheduled for **Phase 3**.

Calling `design.plot_layout()` in Phase 1 raises a `NotImplementedError` explaining that graphical layouts are part of Phase 3, directing users to `design.text_layout()`.

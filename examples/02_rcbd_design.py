"""
Example 02 — RCBD Design Generation
=====================================
Demonstrates creating a Randomized Complete Block Design, printing the block
layout with separators, and exporting to CSV.
"""
import agristatspy as ag

# -------------------------------------------------------------------
# 1. Generate an RCBD with 4 treatments and 3 blocks
# -------------------------------------------------------------------
design = ag.design.rcbd(
    treatments=["Control", "N50", "N100", "N150"],
    replications=3,   # 3 complete blocks
    seed=42,
)

# -------------------------------------------------------------------
# 2. Print text layout with block separators
# -------------------------------------------------------------------
print(design.text_layout(block_separator=True, separator="="))

# -------------------------------------------------------------------
# 3. Validate the design
# -------------------------------------------------------------------
design.validate()
print("Design passed all invariant checks.\n")

# -------------------------------------------------------------------
# 4. Summary
# -------------------------------------------------------------------
print(design.summary())

# -------------------------------------------------------------------
# 5. Export
# -------------------------------------------------------------------
path = design.to_csv("/tmp/rcbd_example.csv")
print(f"\nExported to: {path}")

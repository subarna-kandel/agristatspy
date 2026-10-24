"""
Example 01 — CRD Design Generation
====================================
Demonstrates creating a Completely Randomized Design, inspecting the
assignment table, printing the text layout, and exporting to CSV.
"""
import agristatspy as ag

# -------------------------------------------------------------------
# 1. Generate a CRD with 3 treatments and 4 replications
# -------------------------------------------------------------------
design = ag.design.crd(
    treatments=["Control", "N50", "N100"],
    replications=4,
    seed=12,
)

# -------------------------------------------------------------------
# 2. Inspect the design
# -------------------------------------------------------------------
print("=== Design metadata ===")
print(design.metadata)
print()

print("=== Assignment table ===")
print(design.table.to_string(index=False))
print()

# -------------------------------------------------------------------
# 3. Print text layout (grid mode, default)
# -------------------------------------------------------------------
print("=== Text layout ===")
print(design.text_layout())

# -------------------------------------------------------------------
# 4. Summarize treatment frequencies
# -------------------------------------------------------------------
print("=== Summary ===")
print(design.summary())
print()

# -------------------------------------------------------------------
# 5. Export to CSV
# -------------------------------------------------------------------
path = design.to_csv("/tmp/crd_example.csv")
print(f"Exported to: {path}")

# -------------------------------------------------------------------
# 6. Alternative: auto-generate labels from n_treatments
# -------------------------------------------------------------------
design2 = ag.design.crd(n_treatments=4, replications=5, seed=2026)
print("\nAuto-label design (T1..T4, 5 reps each):")
print(design2.table["treatment"].value_counts().sort_index().to_string())

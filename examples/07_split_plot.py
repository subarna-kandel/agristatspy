"""
Example 07 — Split-plot (irrigation main plots × variety sub-plots)
===================================================================
"""
import agristatspy as ag
from agristatspy.datasets import load_split_plot_example

design = ag.design.split_plot(
    ["Rainfed", "Drip"],
    ["V1", "V2", "V3"],
    replications=3,
    main_name="irrigation",
    sub_name="variety",
    seed=11,
)
print(design.text_layout())
print()

df = load_split_plot_example()
print(ag.anova_split_plot(df, "yield", main="irrigation", sub="variety", block="block"))

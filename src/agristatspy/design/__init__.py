"""
agristatspy.design — experimental-design generation sub-package.

Public interface
----------------
.. code-block:: python

    import agristatspy as ag

    design = ag.design.crd(treatments=["A", "B", "C"], replications=4, seed=1)
    design = ag.design.rcbd(n_treatments=4, replications=3, seed=42)
    design = ag.design.latin_square(n_treatments=4, seed=7)
    design = ag.design.factorial({"N": ["N0", "N50"], "V": ["V1", "V2"]}, replications=3)
    design = ag.design.split_plot(["I0", "I1"], ["V1", "V2", "V3"], replications=3)
    design = ag.design.alpha_lattice(n_treatments=12, block_size=4, replications=3, seed=1)
"""
from agristatspy.design.alpha_lattice import alpha_lattice
from agristatspy.design.crd import crd
from agristatspy.design.factorial import factorial
from agristatspy.design.latin_square import latin_square
from agristatspy.design.rcbd import rcbd
from agristatspy.design.split_plot import split_plot

__all__ = ["crd", "rcbd", "latin_square", "factorial", "split_plot", "alpha_lattice"]

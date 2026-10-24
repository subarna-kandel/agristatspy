"""
agristatspy.design — experimental-design generation sub-package.

Public interface
----------------
.. code-block:: python

    import agristatspy as ag

    design = ag.design.crd(treatments=["A", "B", "C"], replications=4, seed=1)
    design = ag.design.rcbd(n_treatments=4, replications=3, seed=42)
"""
from agristatspy.design.crd import crd
from agristatspy.design.rcbd import rcbd

__all__ = ["crd", "rcbd"]

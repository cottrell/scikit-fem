import numpy as np
from ..element_h1 import ElementH1
from ...refdom import RefTri

class ElementTriCR(ElementH1):
    """The nonconforming Crouzeix-Raviart element."""
    facet_dofs = 1
    maxdeg = 1
    dofnames = ['u']
    doflocs = np.array([[0.5, 0.0], [0.5, 0.5], [0.0, 0.5]])
    refdom = RefTri

    def lbasis(self, X, i):
        (x, y) = X
        if i == 0:
            phi = 1.0 - 2.0 * y
            dphi = np.array([0.0 * x, -2.0 + 0.0 * y])
        elif i == 1:
            phi = 2.0 * x + 2.0 * y - 1.0
            dphi = np.array([2.0 + 0.0 * x, 2.0 + 0.0 * y])
        elif i == 2:
            phi = 1.0 - 2.0 * x
            dphi = np.array([-2.0 + 0.0 * x, 0.0 * x])
        else:
            self._index_error()
        return (phi, dphi)
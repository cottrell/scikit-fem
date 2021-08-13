import numpy as np
from ..element_hdiv import ElementHdiv
from ...refdom import RefTet

class ElementTetRT0(ElementHdiv):
    """The lowest order Raviart-Thomas element."""
    facet_dofs = 1
    maxdeg = 1
    dofnames = ['u^n']
    doflocs = np.array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [0.5, 0.5, 0.5]])
    refdom = RefTet

    def lbasis(self, X, i):
        (x, y, z) = X
        if i == 0:
            phi = np.array([x, y, z - 1.0])
            dphi = 3.0 + 0.0 * x
        elif i == 1:
            phi = np.array([x, y - 1.0, z])
            dphi = 3.0 + 0.0 * x
        elif i == 2:
            phi = np.array([x - 1.0, y, z])
            dphi = 3.0 + 0.0 * x
        elif i == 3:
            phi = np.array([x, y, z])
            dphi = 3.0 + 0.0 * x
        else:
            self._index_error()
        return (phi, dphi)
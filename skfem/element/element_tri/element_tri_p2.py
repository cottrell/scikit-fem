import numpy as np
from ..element_h1 import ElementH1
from ...refdom import RefTri

class ElementTriP2(ElementH1):
    """Piecewise quadratic element."""
    nodal_dofs = 1
    facet_dofs = 1
    maxdeg = 2
    dofnames = ['u', 'u']
    doflocs = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5]])
    refdom = RefTri

    def lbasis(self, X, i):
        (x, y) = X
        if i == 0:
            phi = 1.0 - 3.0 * x - 3.0 * y + 2.0 * x ** 2 + 4.0 * x * y + 2.0 * y ** 2
            dphi = np.array([-3 + 4.0 * x + 4.0 * y, -3 + 4.0 * x + 4.0 * y])
        elif i == 1:
            phi = 2.0 * x ** 2 - x
            dphi = np.array([4.0 * x - 1, 0.0 * x])
        elif i == 2:
            phi = 2.0 * y ** 2 - y
            dphi = np.array([0.0 * x, 4.0 * y - 1])
        elif i == 3:
            phi = 4.0 * x - 4.0 * x ** 2 - 4.0 * x * y
            dphi = np.array([4 - 8.0 * x - 4.0 * y, -4.0 * x])
        elif i == 4:
            phi = 4.0 * x * y
            dphi = np.array([4.0 * y, 4.0 * x])
        elif i == 5:
            phi = 4.0 * y - 4.0 * x * y - 4.0 * y ** 2
            dphi = np.array([-4.0 * y, 4 - 4.0 * x - 8.0 * y])
        else:
            self._index_error()
        return (phi, dphi)
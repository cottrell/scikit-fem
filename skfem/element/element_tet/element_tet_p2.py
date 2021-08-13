import numpy as np
from ..element_h1 import ElementH1
from ...refdom import RefTet

class ElementTetP2(ElementH1):
    """Piecewise quadratic element."""
    nodal_dofs = 1
    edge_dofs = 1
    maxdeg = 2
    dofnames = ['u', 'u']
    doflocs = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    refdom = RefTet

    def lbasis(self, X, i):
        (x, y, z) = X
        if i == 0:
            phi = 1.0 - 3.0 * x + 2.0 * x ** 2 - 3.0 * y + 4.0 * x * y + 2.0 * y ** 2 - 3.0 * z + 4.0 * x * z + 4.0 * y * z + 2.0 * z ** 2
            dphi = np.array([-3.0 + 4.0 * x + 4.0 * y + 4.0 * z, -3.0 + 4.0 * x + 4.0 * y + 4.0 * z, -3.0 + 4.0 * x + 4.0 * y + 4.0 * z])
        elif i == 1:
            phi = -1.0 * x + 2.0 * x ** 2
            dphi = np.array([-1 + 4 * x, 0 * x, 0 * x])
        elif i == 2:
            phi = -1.0 * y + 2.0 * y ** 2
            dphi = np.array([0 * x, -1.0 + 4.0 * y, 0 * x])
        elif i == 3:
            phi = -1.0 * z + 2.0 * z ** 2
            dphi = np.array([0 * x, 0 * x, -1.0 + 4.0 * z])
        elif i == 4:
            phi = 4.0 * x - 4.0 * x ** 2 - 4.0 * x * y - 4 * x * z
            dphi = np.array([4.0 - 8.0 * x - 4.0 * y - 4.0 * z, -4.0 * x, -4.0 * x])
        elif i == 5:
            phi = 4.0 * x * y
            dphi = np.array([4.0 * y, 4.0 * x, 0 * x])
        elif i == 6:
            phi = 0.0 + 4.0 * y - 4.0 * x * y - 4.0 * y ** 2 - 4.0 * y * z
            dphi = np.array([-4.0 * y, 4.0 - 4.0 * x - 8.0 * y - 4.0 * z, -4.0 * y])
        elif i == 7:
            phi = 0.0 + 4.0 * z - 4.0 * x * z - 4.0 * y * z - 4.0 * z ** 2
            dphi = np.array([-4.0 * z, -4.0 * z, 4.0 - 4.0 * x - 4.0 * y - 8.0 * z])
        elif i == 8:
            phi = 0.0 + 4.0 * x * z
            dphi = np.array([4.0 * z, 0 * x, 4 * x])
        elif i == 9:
            phi = 0.0 + 4.0 * y * z
            dphi = np.array([0 * x, 4 * z, 4 * y])
        else:
            self._index_error()
        return (phi, dphi)
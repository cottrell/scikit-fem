from typing import Optional, Tuple
from threading import Thread
from itertools import product
import numpy as np
from numpy import ndarray
from scipy.sparse import csr_matrix
from ..basis import Basis
from .coo_data import COOData
from .form import Form, FormExtraParams

class BilinearForm(Form):
    """A bilinear form for finite element assembly.

    Bilinear forms are defined using functions that takes three arguments:
    trial function ``u``, test function ``v``, and a dictionary of additional
    parameters ``w``.

    >>> from skfem import BilinearForm, InteriorBasis, MeshTri, ElementTriP1
    >>> form = BilinearForm(lambda u, v, _: u * v)
    >>> form.assemble(InteriorBasis(MeshTri(), ElementTriP1())).todense()
    matrix([[0.08333333, 0.04166667, 0.04166667, 0.        ],
            [0.04166667, 0.16666667, 0.08333333, 0.04166667],
            [0.04166667, 0.08333333, 0.16666667, 0.04166667],
            [0.        , 0.04166667, 0.04166667, 0.08333333]])

    Alternatively, you can use :class:`~skfem.assembly.BilinearForm` as a
    decorator:

    >>> @BilinearForm
    ... def form(u, v, _):
    ...     return u * v

    Inside the form definition, ``u`` and ``v`` are tuples containing the basis
    function values at quadrature points.  They also contain the values of
    the derivatives:

    >>> @BilinearForm
    ... def form(u, v, _):
    ...     # u[1][0] is first derivative with respect to x, and so on
    ...     return u[1][0] * v[1][0] + u[1][1] * v[1][1]  # laplacian

    We suggest you to use helper functions from :mod:`skfem.helpers` to make
    the forms more readable:

    >>> from skfem.helpers import dot, grad
    >>> @BilinearForm
    ... def form(u, v, w):
    ...     return dot(grad(u), grad(v))

    """

    def _assemble(self, ubasis: Basis, vbasis: Optional[Basis]=None, **kwargs) -> Tuple[ndarray, ndarray, ndarray, Tuple[int, int]]:
        if vbasis is None:
            vbasis = ubasis
        elif ubasis.X.shape[1] != vbasis.X.shape[1]:
            raise ValueError('Quadrature mismatch: trial and test functions should have same number of integration points.')
        nt = ubasis.nelems
        dx = ubasis.dx
        wdict = FormExtraParams({**ubasis.default_parameters(), **self.dictify(kwargs)})
        sz = ubasis.Nbfun * vbasis.Nbfun * nt
        if self.nthreads > 0:
            data = np.zeros((ubasis.Nbfun, vbasis.Nbfun, nt), dtype=self.dtype)
        else:
            data = np.zeros(sz, dtype=self.dtype)
        rows = np.zeros(sz, dtype=np.int64)
        cols = np.zeros(sz, dtype=np.int64)
        for j in range(ubasis.Nbfun):
            for i in range(vbasis.Nbfun):
                ixs = slice(nt * (vbasis.Nbfun * j + i), nt * (vbasis.Nbfun * j + i + 1))
                rows[ixs] = vbasis.element_dofs[i]
                cols[ixs] = ubasis.element_dofs[j]
                if self.nthreads <= 0:
                    data[ixs] = self._kernel(ubasis.basis[j], vbasis.basis[i], wdict, dx)
        if self.nthreads > 0:
            indices = np.array([[i, j] for (j, i) in product(range(ubasis.Nbfun), range(vbasis.Nbfun))])
            threads = [Thread(target=self._threaded_kernel, args=(data, ix, ubasis.basis, vbasis.basis, wdict, dx)) for ix in np.array_split(indices, self.nthreads, axis=0)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            data = data.flatten('C')
        return (data, rows, cols, (vbasis.N, ubasis.N))

    def coo_data(self, *args, **kwargs) -> COOData:
        return COOData(*self._assemble(*args, **kwargs))

    def assemble(self, *args, **kwargs) -> csr_matrix:
        """Assemble the bilinear form into a sparse matrix.

        Parameters
        ----------
        ubasis
            The :class:`~skfem.assembly.Basis` for ``u``.
        vbasis
            Optionally, specify a different :class:`~skfem.assembly.Basis`
            for ``v``.
        **kwargs
            Any additional keyword arguments are appended to ``w``.

        """
        return COOData._assemble_scipy_csr(*self._assemble(*args, **kwargs))

    def _kernel(self, u, v, w, dx):
        return np.sum(self.form(*u, *v, w) * dx, axis=1)

    def _threaded_kernel(self, data, ix, ubasis, vbasis, wdict, dx):
        for ij in ix:
            (i, j) = ij
            data[j, i] = self._kernel(ubasis[j], vbasis[i], wdict, dx)
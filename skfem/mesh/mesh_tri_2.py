from dataclasses import dataclass, replace
from typing import Type
import numpy as np
from numpy import ndarray
from ..element import Element, ElementTriP2
from .mesh_tri_1 import MeshTri1

@dataclass(repr=False)
class MeshTri2(MeshTri1):
    """A quadratic triangular mesh."""
    doflocs: ndarray = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.0], [0.0, 0.5], [0.5, 0.5], [1.0, 0.5], [0.5, 1.0]], dtype=np.float64).T
    elem: Type[Element] = ElementTriP2
    affine: bool = False
    sort_t: bool = False

    @classmethod
    def init_circle(cls: Type, nrefs: int=3) -> 'MeshTri2':
        m = MeshTri1.init_circle(nrefs=nrefs)
        M = cls.from_mesh(m)
        D = M.dofs.get_facet_dofs(M.boundary_facets()).flatten()
        doflocs = M.doflocs.copy()
        doflocs[:, D] /= np.linalg.norm(doflocs[:, D], axis=0)
        return replace(M, doflocs=doflocs)

    def _repr_svg_(self) -> str:
        from skfem.visuals.svg import draw
        return draw(self, nrefs=2, boundaries_only=True).svg
import warnings
from .abstract_basis import AbstractBasis
from .cell_basis import CellBasis
from .boundary_facet_basis import BoundaryFacetBasis
from .interior_facet_basis import InteriorFacetBasis
from .mortar_facet_basis import MortarFacetBasis
Basis = CellBasis
InteriorBasis = CellBasis
ExteriorFacetBasis = BoundaryFacetBasis

def FacetBasis(*args, side=None, **kwargs):
    """alias of :class:`~skfem.assembly.BoundaryFacetBasis`"""
    if side is None:
        return BoundaryFacetBasis(*args, **kwargs)
    warnings.warn('Initializing FacetBasis using the keyword argument side is deprecated. Use InteriorFacetBasis or MortarFacetBasis instead.', DeprecationWarning)
    if 'mapping' in kwargs:
        if hasattr(kwargs['mapping'], 'helper_to_orig'):
            return MortarFacetBasis(*args, side=side, **kwargs)
    return InteriorFacetBasis(*args, side=side, **kwargs)
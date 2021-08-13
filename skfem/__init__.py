"""Support for wildcard import."""
from .mesh import *
from .assembly import *
from .mapping import *
from .element import *
from .utils import *
from .assembly import __all__ as all_assembly
from .mesh import __all__ as all_mesh
from .element import __all__ as all_element
__all__ = all_mesh + all_assembly + all_element + ['MappingAffine', 'MappingIsoparametric', 'MappingMortar', 'adaptive_theta', 'build_pc_ilu', 'build_pc_diag', 'condense', 'enforce', 'penalize', 'project', 'projection', 'solve', 'solver_direct_scipy', 'solver_eigen_scipy', 'solver_eigen_scipy_sym', 'solver_iter_pcg', 'solver_iter_krylov']
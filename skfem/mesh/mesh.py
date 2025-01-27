from collections import namedtuple
from dataclasses import dataclass, replace
from typing import Callable, Dict, List, Optional, Tuple, Type, Union
from warnings import warn

import numpy as np
from numpy import ndarray

from ..element import BOUNDARY_ELEMENT_MAP, Element


@dataclass(repr=False)
class Mesh:
    """A mesh defined by :class:`~skfem.element.Element` class.

    :class:`~skfem.mesh.Mesh` is defined as a combination of elements/cells by
    specifying the spatial locations of the finite element nodes.

    """

    doflocs: ndarray  #: The locations of the finite element nodes
    t: ndarray  #: The connectivity of the elements/cells
    _boundaries: Optional[Dict[str, ndarray]] = None
    _subdomains: Optional[Dict[str, ndarray]] = None
    elem: Type[Element] = Element
    affine: bool = False
    validate: bool = False  # unused; for backwards compatibility
    # Some parts of the library, most notably the normal vector construction in
    # ElementGlobal._eval_dofs, assume that the element indices are ascending
    # because this leads to consistent normal vectors for both elements sharing
    # a facet.  Therefore, the element indices are sorted in a triangle mesh.
    # However, some algorithms (e.g., adaptive refinement) require switching
    # off this behaviour and, hence, this flag exists.
    sort_t: bool = False

    @property
    def p(self):
        return self.doflocs

    @property
    def dofs(self):
        from skfem.assembly import Dofs
        if not hasattr(self, '_dofs'):
            self._dofs = Dofs(self, self.elem())
        return self._dofs

    @property
    def refdom(self):
        return self.elem.refdom

    @property
    def brefdom(self):
        return self.elem.refdom.brefdom

    @property
    def bndelem(self):
        return BOUNDARY_ELEMENT_MAP[self.elem]()

    @property
    def nelements(self):
        return self.t.shape[1]

    @property
    def nvertices(self):
        return np.max(self.t) + 1

    @property
    def nfacets(self):
        return self.facets.shape[1]

    @property
    def nedges(self):
        return self.edges.shape[1]

    @property
    def nnodes(self):
        return self.t.shape[0]

    @property
    def subdomains(self):
        return self._subdomains

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def facets(self):
        if not hasattr(self, '_facets'):
            self._init_facets()
        return self._facets

    @property
    def t2f(self):
        if not hasattr(self, '_t2f'):
            self._init_facets()
        return self._t2f

    @property
    def f2t(self):
        if not hasattr(self, '_f2t'):
            self._f2t = self.build_inverse(self.t, self.t2f)
        return self._f2t

    @property
    def edges(self):
        if not hasattr(self, '_edges'):
            self._init_edges()
        return self._edges

    @property
    def t2e(self):
        if not hasattr(self, '_t2e'):
            self._init_edges()
        return self._t2e

    def dim(self):
        return self.elem.refdom.dim()

    def boundary_facets(self) -> ndarray:
        """Return an array of boundary facet indices."""
        return np.nonzero(self.f2t[1] == -1)[0]

    def boundary_edges(self) -> ndarray:
        """Return an array of boundary edge indices."""
        facets = self.boundary_facets()
        boundary_edges = np.sort(np.hstack(
            tuple([np.vstack((self.facets[itr, facets],
                              self.facets[(itr + 1) % self.facets.shape[0],
                              facets]))
                   for itr in range(self.facets.shape[0])])).T, axis=1)
        edge_candidates = np.unique(self.t2e[:, self.f2t[0, facets]])
        A = self.edges[:, edge_candidates].T
        B = boundary_edges
        dims = A.max(0) + 1
        ix = np.where(np.in1d(
            np.ravel_multi_index(A.T, dims),  # type: ignore
            np.ravel_multi_index(B.T, dims),  # type: ignore
        ))[0]
        return edge_candidates[ix]

    def define_boundary(self, name: str,
                        test: Callable[[ndarray], ndarray],
                        boundaries_only: bool = True):
        """Define a named boundary via function handle.

        Parameters
        ----------
        name
            Name of the boundary.
        test
            A function which returns True for facet midpoints belonging to the
            boundary.
        boundaries_only
            If True, include only facets on the boundary of the mesh.

        """
        warn("Mesh.define_boundary is deprecated and will be removed in the "
             "next major release.", DeprecationWarning)
        if self._boundaries is None:
            self._boundaries = {}
        self._boundaries[name] = self.facets_satisfying(test, boundaries_only)

    def with_boundaries(self,
                        boundaries: Dict[str, Callable[[ndarray], ndarray]]):
        """Return a copy of the mesh with named boundaries.

        Parameters
        ----------
        boundaries
            A dictionary of lambda functions with the names of the boundaries
            as keys.  The midpoint of the facet should return ``True`` for the
            corresponding lambda function if the facet belongs to the boundary.

        """
        return replace(
            self,
            _boundaries={
                **({} if self._boundaries is None else self._boundaries),
                **{name: self.facets_satisfying(test, True)
                   for name, test in boundaries.items()}
            },
        )

    def with_subdomains(self,
                        subdomains: Dict[str, Callable[[ndarray], ndarray]]):
        """Return a copy of the mesh with named subdomains.

        Parameters
        ----------
        boundaries
            A dictionary of lambda functions with the names of the subdomains
            as keys.  The midpoint of the element should return ``True`` for
            the corresponding lambda function if the element belongs to the
            subdomain.

        """
        return replace(
            self,
            _subdomains={
                **({} if self._subdomains is None else self._subdomains),
                **{name: self.elements_satisfying(test)
                   for name, test in subdomains.items()},
            }
        )

    def boundary_nodes(self) -> ndarray:
        """Return an array of boundary node indices."""
        return np.unique(self.facets[:, self.boundary_facets()])

    def interior_nodes(self) -> ndarray:
        """Return an array of interior node indices."""
        return np.setdiff1d(np.arange(0, self.p.shape[1]),
                            self.boundary_nodes())

    def nodes_satisfying(self,
                         test: Callable[[ndarray], ndarray],
                         boundaries_only: bool = False) -> ndarray:
        """Return nodes that satisfy some condition.

        Parameters
        ----------
        test
            A function which returns ``True`` for the set of nodes that are to
            be included in the return set.
        boundaries_only
            If ``True``, include only boundary facets.

        """
        nodes = np.nonzero(test(self.p))[0]
        if boundaries_only:
            nodes = np.intersect1d(nodes, self.boundary_nodes())
        return nodes

    def facets_satisfying(self,
                          test: Callable[[ndarray], ndarray],
                          boundaries_only: bool = False) -> ndarray:
        """Return facets whose midpoints satisfy some condition.

        Parameters
        ----------
        test
            A function which returns ``True`` for the facet midpoints that are
            to be included in the return set.
        boundaries_only
            If ``True``, include only boundary facets.

        """
        midp = [np.sum(self.p[itr, self.facets], axis=0) / self.facets.shape[0]
                for itr in range(self.dim())]
        facets = np.nonzero(test(np.array(midp)))[0]
        if boundaries_only:
            facets = np.intersect1d(facets, self.boundary_facets())
        return facets

    def elements_satisfying(self,
                            test: Callable[[ndarray], ndarray]) -> ndarray:
        """Return elements whose midpoints satisfy some condition.

        Parameters
        ----------
        test
            A function which returns ``True`` for the element midpoints that
            are to be included in the return set.

        """
        midp = [np.sum(self.p[itr, self.t], axis=0) / self.t.shape[0]
                for itr in range(self.dim())]
        return np.nonzero(test(np.array(midp)))[0]

    def _expand_facets(self, ix: ndarray) -> Tuple[ndarray, ndarray]:
        """Return vertices and edges corresponding to given facet indices.

        Parameters
        ----------
        ix
            An array of facet indices.

        """
        vertices = np.unique(self.facets[:, ix].flatten())

        if self.dim() == 3:
            edge_candidates = self.t2e[:, self.f2t[0, ix]].flatten()
            # subset of edges that share all points with the given facets
            subset = np.nonzero(
                np.prod(np.isin(self.edges[:, edge_candidates],
                                self.facets[:, ix].flatten()),
                        axis=0)
            )[0]
            edges = np.intersect1d(self.boundary_edges(),
                                   edge_candidates[subset])
        else:
            edges = np.array([], dtype=np.int64)

        return vertices, edges

    def _mapping(self):
        """Return a default reference mapping for the mesh."""
        from skfem.mapping import MappingAffine, MappingIsoparametric
        if not hasattr(self, '_cached_mapping'):
            if self.affine:
                self._cached_mapping = MappingAffine(self)
            else:
                # TODO make MappingIsoparametric compatible with self
                FakeMesh = namedtuple(
                    'FakeMesh',
                    ['p', 't', 'facets', 't2f', 'f2t', 'dim']
                )
                fakemesh = FakeMesh(
                    self.doflocs,
                    self.dofs.element_dofs,
                    self.facets,
                    self.t2f,
                    self.f2t,
                    lambda: self.dim(),
                )
                self._cached_mapping = MappingIsoparametric(
                    fakemesh,
                    self.elem(),
                    self.bndelem,
                )
        return self._cached_mapping

    def _init_facets(self):
        """Initialize ``self.facets``."""
        self._facets, self._t2f = self.build_entities(
            self.t,
            self.elem.refdom.facets,
        )

    def _init_edges(self):
        """Initialize ``self.edges``."""
        self._edges, self._t2e = self.build_entities(
            self.t,
            self.elem.refdom.edges,
        )

    def __post_init__(self):
        """Support node orders used in external formats.

        We expect ``self.doflocs`` to be ordered based on the
        degrees-of-freedom in :class:`skfem.assembly.Dofs`.  External formats
        for high order meshes commonly use a less strict ordering scheme and
        the extra nodes are described as additional rows in ``self.t``.  This
        method attempts to accommodate external formas by reordering
        ``self.doflocs`` and changing the indices in ``self.t``.

        """
        if self.sort_t:
            self.t = np.sort(self.t, axis=0)

        self.doflocs = np.asarray(self.doflocs, dtype=np.float64, order="K")
        self.t = np.asarray(self.t, dtype=np.int64, order="K")

        M = self.elem.refdom.nnodes

        if self.nnodes > M:
            # reorder DOFs to the expected format: vertex DOFs are first
            p, t = self.doflocs, self.t
            t_nodes = t[:M]
            uniq, ix = np.unique(t_nodes, return_inverse=True)
            self.t = (np.arange(len(uniq), dtype=np.int64)[ix]
                      .reshape(t_nodes.shape))
            doflocs = np.hstack((
                p[:, uniq],
                np.zeros((p.shape[0], np.max(t) + 1 - len(uniq))),
            ))
            doflocs[:, self.dofs.element_dofs[M:].flatten('F')] =\
                p[:, t[M:].flatten('F')]
            self.doflocs = doflocs

        # C_CONTIGUOUS is more performant in dimension-based slices
        if not self.doflocs.flags['C_CONTIGUOUS']:
            if self.doflocs.shape[1] > 1000:
                warn("Transforming over 1000 vertices to C_CONTIGUOUS.")
            self.doflocs = np.ascontiguousarray(self.doflocs)

        if not self.t.flags['C_CONTIGUOUS']:
            if self.t.shape[1] > 1000:
                warn("Transforming over 1000 elements to C_CONTIGUOUS.")
            self.t = np.ascontiguousarray(self.t)

    def __add__(self, other):
        """Join two meshes."""
        if not isinstance(other, type(self)):
            raise TypeError("Can only join meshes with same type.")
        p = np.hstack((self.p, other.p))
        t = np.hstack((self.t, other.t + self.p.shape[1]))
        tmp = np.ascontiguousarray(p.T)
        tmp, ixa, ixb = np.unique(tmp.view([('', tmp.dtype)] * tmp.shape[1]),
                                  return_index=True, return_inverse=True)
        p = p[:, ixa]
        t = ixb[t]
        cls = type(self)
        return cls(p, t)

    def __repr__(self):
        return "{} mesh with {} vertices and {} elements.".format(
            self.elem.refdom.name,
            self.nvertices,
            self.nelements,
        )

    def __str__(self):
        return self.__repr__()

    def save(self,
             filename: str,
             point_data: Optional[Dict[str, ndarray]] = None,
             **kwargs) -> None:
        """Export the mesh and fields using meshio.

        Parameters
        ----------
        filename
            The output filename, with suffix determining format;
            e.g. .msh, .vtk, .xdmf
        point_data
            Data related to the vertices of the mesh.

        """
        from skfem.io.meshio import to_file
        return to_file(self, filename, point_data, **kwargs)

    @classmethod
    def load(cls, filename: str):
        """Load a mesh using meshio.

        Parameters
        ----------
        filename
            The filename of the mesh file.

        """
        from skfem.io.meshio import from_file
        return from_file(filename)

    @classmethod
    def from_dict(cls, data):
        """For backwards compatibility."""
        if 'p' not in data or 't' not in data:
            raise ValueError("Dictionary must contain keys 'p' and 't'.")
        else:
            data['p'] = np.ascontiguousarray(np.array(data['p']).T)
            data['t'] = np.ascontiguousarray(np.array(data['t']).T)
        if 'boundaries' in data and data['boundaries'] is not None:
            data['boundaries'] = {k: np.array(v)
                                  for k, v in data['boundaries'].items()}
        if 'subdomains' in data and data['subdomains'] is not None:
            data['subdomains'] = {k: np.array(v)
                                  for k, v in data['subdomains'].items()}
        data['doflocs'] = data.pop('p')
        data['_subdomains'] = data.pop('subdomains')
        data['_boundaries'] = data.pop('boundaries')
        return cls(**data)

    def to_dict(self) -> Dict[str, Optional[Dict[str, List[float]]]]:
        """For backwards compatibility."""
        boundaries: Optional[Dict[str, List[float]]] = None
        subdomains: Optional[Dict[str, List[float]]] = None
        if self.boundaries is not None:
            boundaries = {k: v.tolist() for k, v in self.boundaries.items()}
        if self.subdomains is not None:
            subdomains = {k: v.tolist() for k, v in self.subdomains.items()}
        return {
            'p': self.p.T.tolist(),
            't': self.t.T.tolist(),
            'boundaries': boundaries,
            'subdomains': subdomains,
        }

    @classmethod
    def from_mesh(cls, mesh):
        """Reuse an existing mesh by adding nodes.

        Parameters
        ----------
        mesh
            The mesh used in the initialization.  Connectivity of the new mesh
            will match ``mesh.t``.

        """
        from skfem.assembly import Dofs

        mapping = mesh._mapping()
        nelem = cls.elem
        dofs = Dofs(mesh, nelem())
        locs = mapping.F(nelem.doflocs.T)
        doflocs = np.zeros((locs.shape[0], dofs.N))

        # match mapped dofs and global dof numbering
        for itr in range(locs.shape[0]):
            for jtr in range(dofs.element_dofs.shape[0]):
                doflocs[itr, dofs.element_dofs[jtr]] = locs[itr, :, jtr]

        return cls(
            doflocs=doflocs,
            t=mesh.t,
        )

    @classmethod
    def init_refdom(cls):
        """Initialize a mesh corresponding to the reference domain."""
        return cls(cls.elem.refdom.p, cls.elem.refdom.t)

    def refined(self, times_or_ix: Union[int, ndarray] = 1):
        """Return a refined mesh.

        Parameters
        ----------
        times_or_ix
            Either an integer giving the number of uniform refinements or an
            array of element indices for adaptive refinement.

        """
        m = self
        if isinstance(times_or_ix, int):
            for _ in range(times_or_ix):
                m = m._uniform()
        else:
            m = m._adaptive(times_or_ix)
        return m

    def scaled(self, factors):
        """Return a new mesh with scaled dimensions.

        Parameters
        ----------
        factors
            Scale each dimension by a factor.

        """
        if isinstance(factors, float):
            # for backwards compatibility
            factors = self.doflocs.shape[0] * [factors]
        return replace(
            self,
            doflocs=np.array([self.doflocs[itr] * factors[itr]
                              for itr in range(len(factors))]),
        )

    def translated(self, diffs):
        """Return a new translated mesh.

        Parameters
        ----------
        diffs
            Translate the mesh by a vector. Must have same size as the mesh
            dimension.

        """
        return replace(
            self,
            doflocs=np.array([self.doflocs[itr] + diffs[itr]
                              for itr in range(len(diffs))]),
        )

    def mirrored(self,
                 normal: Tuple[float, ...],
                 point: Optional[Tuple[float, ...]] = None):
        """Return a mesh mirrored with respect to a normal.

        Meant to be combined with the other methods to build more general
        meshes, e.g.,

        >>> from skfem import MeshTet
        >>> m1 = MeshTet()
        >>> m2 = m1.mirrored((1, 0, 0))
        >>> m3 = m1.mirrored((0, 1, 0))
        >>> m4 = m1.mirrored((0, 0, 1))
        >>> m = m1 + m2 + m3 + m4
        >>> (m.nvertices, m.nelements)
        (20, 20)

        Parameters
        ----------
        normal
            The normal vector of the mirror plane.
        point
            An optional point through which the plane passes. By default, the
            point corresponds to the origin.

        """
        if point is None:
            point = (0,) * self.dim()

        p = self.p.copy()
        p0 = np.array(point)
        n = np.array(normal)
        n = n / np.linalg.norm(n)
        p = p - 2. * np.dot(n, p - p0[:, None]) * n[:, None]

        return replace(
            self,
            doflocs=p,
        )

    def _uniform(self):
        """Perform a single uniform refinement."""
        raise NotImplementedError

    def _adaptive(self, ix: ndarray):
        """Adaptively refine the given set of elements."""
        raise NotImplementedError

    def _splitref(self, nrefs: int = 1):
        """Split mesh into separate nonconnected elements and refine.

        Used for visualization purposes.

        Parameters
        ----------
        nrefs
            The number of refinements.

        """
        cls = type(self)
        m = cls.init_refdom().refined(nrefs)
        X = m.p
        x = self._mapping().F(m.p)

        # create connectivity for the new mesh
        nt = self.nelements
        t = np.tile(m.t, (1, nt))
        dt = np.max(t)
        t += ((dt + 1)
              * (np.tile(np.arange(nt), (m.t.shape[0] * m.t.shape[1], 1))
                 .flatten('F')
                 .reshape((-1, m.t.shape[0])).T))

        if X.shape[0] == 1:
            p = np.array([x.flatten()])
        else:
            p = x[0].flatten()
            for itr in range(len(x) - 1):
                p = np.vstack((p, x[itr + 1].flatten()))

        return cls(p, t)

    @staticmethod
    def build_entities(t, indices, sort=True):
        """Build low dimensional topological entities."""
        indexing = np.hstack(tuple([t[ix] for ix in indices]))
        sorted_indexing = np.sort(indexing, axis=0)

        sorted_indexing, ixa, ixb = np.unique(sorted_indexing,
                                              axis=1,
                                              return_index=True,
                                              return_inverse=True)
        mapping = ixb.reshape((len(indices), t.shape[1]))

        if sort:
            return np.ascontiguousarray(sorted_indexing), mapping

        return np.ascontiguousarray(indexing[:, ixa]), mapping

    @staticmethod
    def build_inverse(t, mapping):
        """Build inverse mapping from low dimensional topological entities."""
        e = mapping.flatten(order='C')
        tix = np.tile(np.arange(t.shape[1]), (1, t.shape[0]))[0]

        e_first, ix_first = np.unique(e, return_index=True)
        e_last, ix_last = np.unique(e[::-1], return_index=True)
        ix_last = e.shape[0] - ix_last - 1

        inverse = np.zeros((2, np.max(mapping) + 1), dtype=np.int64)
        inverse[0, e_first] = tix[ix_first]
        inverse[1, e_last] = tix[ix_last]
        inverse[1, np.nonzero(inverse[0] == inverse[1])[0]] = -1

        return inverse

    @staticmethod
    def strip_extra_coordinates(p: ndarray) -> ndarray:
        """Fallback for 3D meshes."""
        return p

    def param(self) -> float:
        """Return mesh parameter, viz the length of the longest edge."""
        raise NotImplementedError

    def _reix(self, ix: ndarray) -> Tuple[ndarray, ndarray]:
        """Connect ``self.p`` based on the indices ``ix``."""
        ixuniq = np.unique(ix)
        t = np.zeros(np.max(ix) + 1, dtype=np.int64)
        t[ixuniq] = np.arange(len(ixuniq), dtype=np.int64)
        return self.p[:, ixuniq], t[ix]

    def remove_elements(self, element_indices: ndarray):
        """Construct a new mesh by removing elements.

        Parameters
        ----------
        element_indices
            List of element indices to remove.

        """
        p, t = self._reix(np.delete(self.t, element_indices, axis=1))
        return replace(
            self,
            doflocs=p,
            t=t,
        )

    def element_finder(self, mapping=None):
        """Return a function handle from location to element index.

        Parameters
        ----------
        mapping
            The affine mapping for the mesh.

        """
        raise NotImplementedError

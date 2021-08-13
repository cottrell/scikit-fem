"""Import/export any formats supported by meshio."""
import meshio
import numpy as np
from ... import skfem
MESH_TYPE_MAPPING = {'tetra': skfem.MeshTet, 'hexahedron': skfem.MeshHex, 'triangle': skfem.MeshTri, 'quad': skfem.MeshQuad, 'line': skfem.MeshLine, 'tetra10': skfem.MeshTet2, 'triangle6': skfem.MeshTri2, 'quad9': skfem.MeshQuad2}
BOUNDARY_TYPE_MAPPING = {'line': 'vertex', 'triangle': 'line', 'quad': 'line', 'tetra': 'triangle', 'hexahedron': 'quad', 'tetra10': 'triangle', 'triangle6': 'line', 'quad9': 'line'}
TYPE_MESH_MAPPING = {MESH_TYPE_MAPPING[k]: k for k in dict(reversed(list(MESH_TYPE_MAPPING.items())))}

def from_meshio(m, out=None, int_data_to_sets=False):
    cells = m.cells_dict
    meshio_type = None
    for k in cells:
        if k in {'tetra', 'hexahedron', 'tetra10'}:
            meshio_type = k
            break
    if meshio_type is None:
        for k in cells:
            if k in {'triangle', 'quad', 'triangle6', 'quad9'}:
                meshio_type = k
                break
    if meshio_type is None:
        for k in cells:
            if k == 'line':
                meshio_type = k
                break
    if meshio_type is None:
        raise NotImplementedError('Mesh type(s) not supported in import: {}.'.format(cells.keys()))
    mesh_type = MESH_TYPE_MAPPING[meshio_type]
    p = np.ascontiguousarray(mesh_type.strip_extra_coordinates(m.points).T)
    t = np.ascontiguousarray(cells[meshio_type].T)
    if meshio_type == 'hexahedron':
        t = t[[0, 4, 3, 1, 7, 5, 2, 6]]
    if int_data_to_sets:
        m.int_data_to_sets()
    subdomains = {}
    boundaries = {}
    if m.cell_sets:
        subdomains = {k: v[meshio_type] for (k, v) in m.cell_sets_dict.items() if meshio_type in v}
    mtmp = mesh_type(p, t)
    bnd_type = BOUNDARY_TYPE_MAPPING[meshio_type]
    if m.cell_sets and bnd_type in m.cells_dict:
        facets = {k: [tuple(f) for f in np.sort(m.cells_dict[bnd_type][v[bnd_type]])] for (k, v) in m.cell_sets_dict.items() if bnd_type in v}
        boundaries = {k: np.array([i for (i, f) in enumerate(map(tuple, mtmp.facets.T)) if f in v]) for (k, v) in facets.items()}
    if m.cell_data and m.field_data:
        try:
            elements_tag = m.cell_data_dict['gmsh:physical'][meshio_type]
            subdomains = {}
            tags = np.unique(elements_tag)

            def find_tagname(tag):
                for key in m.field_data:
                    if m.field_data[key][0] == tag:
                        return key
                return None
            for tag in tags:
                t_set = np.nonzero(tag == elements_tag)[0]
                subdomains[find_tagname(tag)] = t_set
            if bnd_type in m.cell_data_dict['gmsh:physical']:
                facets = m.cells_dict[bnd_type]
                facets_tag = m.cell_data_dict['gmsh:physical'][bnd_type]
            dic = {tuple(np.sort(facets[i])): facets_tag[i] for i in range(facets.shape[0])}
            index = np.array([[dic[tuple(np.sort(mtmp.facets[:, i]))], i] for i in mtmp.boundary_facets() if tuple(np.sort(mtmp.facets[:, i])) in dic])
            tags = index[:, 0]
            boundaries = {}
            for tag in np.unique(tags):
                tagindex = np.nonzero(tags == tag)[0]
                boundaries[find_tagname(tag)] = index[tagindex, 1]
        except Exception:
            pass
    if m.cell_data:
        (_boundaries, _subdomains) = mtmp._decode_cell_data(m.cell_data)
        boundaries.update(_boundaries)
        subdomains.update(_subdomains)
    if out is not None and isinstance(out, list):
        for (i, field) in enumerate(out):
            out[i] = getattr(m, field)
    return mesh_type(p, t, None if len(boundaries) == 0 else boundaries, None if len(subdomains) == 0 else subdomains)

def from_file(filename, out, **kwargs):
    return from_meshio(meshio.read(filename), out, **kwargs)

def to_meshio(mesh, point_data=None, cell_data=None):
    t = mesh.dofs.element_dofs.copy()
    if isinstance(mesh, skfem.MeshHex):
        t = t[[0, 3, 6, 2, 1, 5, 7, 4]]
    mtype = TYPE_MESH_MAPPING[type(mesh)]
    cells = {mtype: t.T}
    if cell_data is None:
        cell_data = {}
    cell_data.update(mesh._encode_cell_data())
    mio = meshio.Mesh(mesh.p.T, cells, point_data=point_data, cell_data=cell_data)
    return mio

def to_file(mesh, filename, point_data=None, cell_data=None, **kwargs):
    meshio.write(filename, to_meshio(mesh, point_data, cell_data), **kwargs)
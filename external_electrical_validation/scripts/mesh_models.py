"""Body-fitted tetrahedral meshes. Coordinates are mm; conversion is in driver."""
from pathlib import Path
import json
import numpy as np
import gmsh


def _collect(fracture_tags=()):
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(tags)
    ids = np.asarray(tags)[order]
    nodes = np.asarray(xyz).reshape(-1, 3)[order]
    cells, labels = [], []
    for _, tag in gmsh.model.getEntities(3):
        types, _, ntags = gmsh.model.mesh.getElements(3, tag)
        for typ, ns in zip(types, ntags):
            if typ != 4:
                raise RuntimeError(f"Expected linear tetrahedra, got {typ}")
            c = np.searchsorted(ids, np.asarray(ns).reshape(-1, 4))
            cells.append(c)
            labels.extend([int(tag in fracture_tags)] * len(c))
    return nodes, np.vstack(cells), np.asarray(labels, dtype=np.int8)


def cell_mesh(path, alpha=0.021, spacing=3.0, level=1, epsilon=0.062):
    """One oblate spheroid in a volume fixed by measured crack density.

    This is an explicitly assumed regular cell, not reconstructed crack positions.
    Opposite principal-axis faces are electrodes; transverse faces are insulated.
    """
    path = Path(path)
    if path.exists():
        return dict(np.load(path))
    a = 1.5
    volume = a**3 / epsilon
    width = np.sqrt(volume / spacing)
    lengths = np.array([width, width, spacing])
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("General.NumThreads", 6)
        gmsh.model.add("single_penny_cell")
        # Mesh a round sphere in stretched coordinates, then map all z nodes
        # linearly by alpha. This gives the exact oblate geometry and avoids
        # singular surface parametrization of an OCC-dilated B-spline sphere.
        computational_lengths = lengths / np.array([1, 1, alpha])
        box = gmsh.model.occ.addBox(*(-computational_lengths/2), *computational_lengths)
        ball = gmsh.model.occ.addSphere(0, 0, 0, a)
        _, maps = gmsh.model.occ.fragment([(3, box)], [(3, ball)])
        fracture = [tag for dim, tag in maps[1] if dim == 3]
        gmsh.model.occ.synchronize()
        print('OCC cell geometry ready', flush=True)
        hs = [0.24, 0.16, 0.105, 0.07][level]
        hf = [1.5, 1.0, 0.7, 0.48][level]
        surf = sorted(set(tag for d, tag in gmsh.model.getBoundary([(3,t) for t in fracture], oriented=False) if d == 2))
        dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist,"SurfacesList",surf)
        gmsh.model.mesh.field.setNumber(dist,"Sampling",80)
        th = gmsh.model.mesh.field.add("Threshold")
        for key, val in {"InField":dist,"SizeMin":hs,"SizeMax":hf,"DistMin":0.08,"DistMax":0.7}.items():
            gmsh.model.mesh.field.setNumber(th,key,val)
        gmsh.model.mesh.field.setAsBackgroundMesh(th)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", hs*0.20)
        gmsh.option.setNumber("Mesh.MeshSizeMax", hf)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize()
        nodes, cells, labels = _collect(fracture)
        nodes[:, 2] *= alpha
        result = dict(nodes=nodes, cells=cells, labels=labels, lengths=lengths,
                      geometric_volume=np.array(volume), alpha=np.array(alpha))
        path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path, **result)
        return result
    finally:
        gmsh.finalize()


def cylinder_mesh(path, length=40., diameter=25.4, level=1):
    path = Path(path)
    if path.exists():
        return dict(np.load(path))
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 6)
        gmsh.model.add("end_electrode_specimen")
        gmsh.model.occ.addCylinder(0,0,0,0,0,length,diameter/2)
        gmsh.model.occ.synchronize()
        h = [3.,2.,1.3,0.9][level]
        gmsh.option.setNumber("Mesh.MeshSizeMin",h)
        gmsh.option.setNumber("Mesh.MeshSizeMax",h)
        gmsh.model.mesh.generate(3)
        nodes,cells,labels = _collect()
        result=dict(nodes=nodes,cells=cells,labels=labels,length=np.array(length),diameter=np.array(diameter))
        path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,**result)
        return result
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    import time
    root=Path(__file__).resolve().parents[1]
    start=time.time()
    m=cell_mesh(root/"runs"/"meshes"/"cell_a0021_s3_l0.npz",level=0)
    xyz=m['nodes'][m['cells']]
    vols=np.abs(np.linalg.det(xyz[:,1:]-xyz[:,:1]))/6
    print(json.dumps(dict(nodes=len(m['nodes']),tets=len(m['cells']),volume=vols.sum(),
                         crack_fraction=vols[m['labels']==1].sum()/vols.sum(),seconds=time.time()-start)))

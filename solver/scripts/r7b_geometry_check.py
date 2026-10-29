#!/usr/bin/env python3
"""Read-only measurements of an actual solver ellipsoid boundary, in metres.

This module NEVER builds a fracture or relabels a tetrahedron.  The caller must
export the real material-union boundary, not an analytic display surface.
See docs/r7b_geometry_check_interface_CN.md for MATLAB export and input schema.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components


U = np.array([.5, 0., math.sqrt(3.) / 2.])
V = np.array([0., 1., 0.])
N = np.array([-math.sqrt(3.) / 2., 0., .5])
BASIS = np.column_stack([U, V, N])
SPACING = {"S1": 8.683e-6, "S5": 3.997e-6}


def definition(sample: str, parameters: dict | None = None) -> dict[str, Any]:
    """Actual declared scale, with an explicit technical/nonformal marker.

    No inference from mesh extents and no automatic rescaling of measured
    points.  Parameter-file dimensions are independently checked for internal
    consistency and the prescribed unchanged ellipsoid proportions/rotation.
    """
    p = parameters or {}
    if parameters is None:
        h, L = SPACING[sample], 128 * SPACING[sample]
    else:
        if "sample" in p and str(p["sample"]) != sample:
            raise ValueError("Parameter sample conflicts with --sample")
        h = p.get("spacing_m", p.get("voxel_size_m"))
        L = p.get("L_m", p.get("sample_L_m"))
        if h is None or L is None:
            raise ValueError("Explicit parameters require spacing_m/voxel_size_m AND L_m/sample_L_m")
        h, L = float(h), float(L)
        if not np.isfinite(h) or not np.isfinite(L) or h <= 0 or L <= 0:
            raise ValueError("Positive finite actual metre spacing and sample length required")
        if not np.isclose(h, SPACING[sample], rtol=1e-12, atol=0):
            raise ValueError("Actual CT voxel spacing must remain the original sample spacing")
    shape = np.asarray(p.get("shape_xyz", [round(L/h)]*3), dtype=int).ravel()
    if shape.shape != (3,) or not np.all(shape == shape[0]) or not np.isclose(L/h, shape[0], rtol=1e-12, atol=1e-10):
        raise ValueError("Actual cubic shape, spacing and L are inconsistent")
    axes = np.asarray(p.get("ellipsoid_semiaxes_m", p.get("semiaxes_m", L*np.array([.25, .25, .0125]))), dtype=float).ravel()
    center = np.asarray(p.get("ellipsoid_center_m", p.get("center_m", np.full(3, L/2))), dtype=float).ravel()
    if axes.shape != (3,) or not np.allclose(axes/L, [.25, .25, .0125], rtol=1e-12, atol=1e-14):
        raise ValueError("Parameters changed the prescribed ellipsoid half-axis ratios")
    if center.shape != (3,) or not np.allclose(center/L, [.5, .5, .5], rtol=1e-12, atol=1e-14):
        raise ValueError("Parameters changed the prescribed ellipsoid centre")
    for key, expected in (("basis_u", U), ("basis_v", V), ("basis_n", N)):
        if key in p and not np.allclose(np.asarray(p[key]).ravel(), expected, rtol=1e-12, atol=1e-14):
            raise ValueError("Parameters changed the prescribed 60-degree orientation")
    technical = bool(shape[0] != 128 or p.get("technical_pilot_only", False) or p.get("mode") == "technical32")
    scope = f"technical_native_centre{shape[0]}_NOT_FORMAL_CT" if technical else "prescribed_native_centre128_geometry"
    if p.get("mode") == "uniform_pilot":
        technical = True
        scope = "uniform_native_ellipsoid_INTERFACE_PILOT_NOT_CT"
    radius = np.sqrt((BASIS**2) @ (axes**2))
    return dict(sample=sample, L_m=L, voxel_size_m=h, shape_xyz=shape,
                geometry_scope=scope, technical_pilot_only=technical,
                formal_CT_scientific_result=False,
                source_parameters_used=parameters is not None,
                protection_regions_present=bool(p.get("protection_regions_present", not technical)),
                center_m=center, semiaxes_m=axes, basis=BASIS,
                ideal_volume_m3=4 * np.pi * np.prod(axes) / 3,
                ideal_bbox_min_m=center-radius, ideal_bbox_max_m=center+radius)


def _rows(array: np.ndarray, width: int, name: str) -> np.ndarray:
    a = np.asarray(array)
    if a.ndim != 2:
        raise ValueError(f"{name} must be a 2-D array")
    if a.shape[1] != width and a.shape[0] == width:
        a = a.T
    if a.shape[1] != width:
        raise ValueError(f"{name} needs {width} columns; got {a.shape}")
    return a


def _dump(path: Path, data: Any) -> None:
    def cast(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, default=cast,
                               allow_nan=False) + "\n", encoding="utf-8")


def _csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def spheroid_signed_distance(points_local_m: np.ndarray, a: float, c: float) -> np.ndarray:
    """Signed nearest Euclidean distance to a=b spheroid (negative inside).

    The nearest-point displacement is normal to the analytic surface.  Solve
    its secular equation, including the equatorial medial-axis special case.
    This measures actual face points; it does not project/update the mesh.
    """
    q = np.asarray(points_local_m, dtype=float)
    radial = np.hypot(q[:, 0], q[:, 1])
    z = np.abs(q[:, 2])
    implicit = (radial / a)**2 + (z / c)**2 - 1
    closest_r = np.empty_like(radial)
    closest_z = np.empty_like(z)
    equator = z <= a * 1e-14
    medial = equator & (radial < (a*a-c*c)/a)
    closest_r[medial] = a*a * radial[medial] / (a*a-c*c)
    closest_z[medial] = c * np.sqrt(np.maximum(0, 1-(closest_r[medial]/a)**2))
    regular_eq = equator & ~medial
    closest_r[regular_eq] = a
    closest_z[regular_eq] = 0
    take = ~equator
    r = radial[take] / c
    zz = z[take] / c
    d = (a/c)**2
    lo = np.full(r.shape, -1.)
    hi = np.maximum(np.sqrt(d)*r + zz, 1.)
    for _ in range(68):
        t = (lo + hi) / 2
        f = d*r*r / (d+t)**2 + zz*zz / (1+t)**2 - 1
        lo = np.where(f > 0, t, lo)
        hi = np.where(f > 0, hi, t)
    t = (lo+hi)/2
    closest_r[take] = c*d*r/(d+t)
    closest_z[take] = c*zz/(1+t)
    dist = np.hypot(closest_r-radial, closest_z-z)
    return np.where(implicit >= 0, dist, -dist)


def face_topology(faces: np.ndarray) -> dict:
    f = np.asarray(faces, dtype=np.int64)
    # An interface shared by several CT subdomains must occur only ONCE here.
    edges = np.sort(np.concatenate((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]])), axis=1)
    unique_edges, inverse, counts = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind="stable")
    starts = np.r_[0, np.cumsum(counts[:-1])]
    manifold = counts == 2
    edge_faces = np.tile(np.arange(len(f), dtype=np.int64), 3)
    left = edge_faces[order[starts[manifold]]]
    right = edge_faces[order[starts[manifold]+1]]
    graph = sparse.coo_matrix((np.ones(len(left), dtype=np.uint8), (left, right)),
                              shape=(len(f), len(f))).tocsr()
    component_count, labels = connected_components(graph, directed=False)
    component_sizes = np.bincount(labels)
    unique_faces = np.unique(np.sort(f, axis=1), axis=0)
    result = dict(face_count=len(f), used_vertex_count=len(np.unique(f)),
                  unique_edge_count=len(unique_edges),
                  surface_edge_connected_components=int(component_count),
                  surface_component_triangle_counts=sorted(component_sizes.tolist(), reverse=True),
                  open_edges=int(np.count_nonzero(counts == 1)),
                  nonmanifold_edges=int(np.count_nonzero(counts > 2)),
                  duplicate_faces=int(len(f)-len(unique_faces)),
                  euler_characteristic=int(len(np.unique(f))-len(unique_edges)+len(f)),
                  closed_two_manifold=bool(np.all(counts == 2) and len(unique_faces) == len(f)))
    return result


def tetrahedron_checks(vertices_m: np.ndarray, tets_inside: np.ndarray) -> tuple[dict, np.ndarray]:
    """Only actual, already fitted inside tets; never determines inside labels."""
    t = np.asarray(tets_inside, dtype=np.int64)
    volume = 0.
    origin = vertices_m[t[0, 0]].copy()
    first_moment = np.zeros(3)
    second_moment = np.zeros((3, 3))
    zero_count = 0
    signed_negative = 0
    for start in range(0, len(t), 100_000):
        p = vertices_m[t[start:start+100_000]]
        det = np.einsum("ij,ij->i", p[:, 1]-p[:, 0],
                        np.cross(p[:, 2]-p[:, 0], p[:, 3]-p[:, 0]))
        volume += np.sum(np.abs(det))/6
        vv = np.abs(det)/6
        rel = p-origin
        sums = rel.sum(axis=1)
        first_moment += np.einsum("i,ij->j", vv/4, sums)
        second_moment += (np.einsum("i,ij,ik->jk", vv/20, sums, sums)
                          + np.einsum("i,ilj,ilk->jk", vv/20, rel, rel))
        zero_count += int(np.count_nonzero(det == 0))
        signed_negative += int(np.count_nonzero(det < 0))
    all_faces = np.sort(np.concatenate((t[:, [1, 2, 3]], t[:, [0, 2, 3]],
                                       t[:, [0, 1, 3]], t[:, [0, 1, 2]])), axis=1)
    unique, inv, count = np.unique(all_faces, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inv, kind="stable")
    starts = np.r_[0, np.cumsum(count[:-1])]
    owner = np.tile(np.arange(len(t), dtype=np.int64), 4)
    take = count == 2
    left = owner[order[starts[take]]]
    right = owner[order[starts[take]+1]]
    graph = sparse.coo_matrix((np.ones(len(left), dtype=np.uint8), (left, right)),
                              shape=(len(t), len(t))).tocsr()
    nc, label = connected_components(graph, directed=False)
    centroid_rel = first_moment/volume
    covariance = second_moment/volume-np.outer(centroid_rel, centroid_rel)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    equivalent_axes = np.sqrt(np.maximum(5*eigenvalues, 0))
    thin_direction = eigenvectors[:, 0]
    if np.dot(thin_direction, N) < 0:
        thin_direction = -thin_direction
    thin_angle = np.degrees(np.arccos(np.clip(np.dot(thin_direction, N), -1, 1)))
    return dict(inside_tetrahedron_count=len(t), actual_tetrahedron_volume_m3=float(volume),
                actual_uniform_volume_centroid_m=origin+centroid_rel,
                actual_uniform_volume_covariance_m2=covariance,
                moment_equivalent_semiaxes_ascending_m=equivalent_axes,
                moment_equivalent_thin_axis_unit=thin_direction,
                moment_thin_axis_angle_from_prescribed_n_deg=float(thin_angle),
                moment_axes_are_not_exact_boundary_extents=True,
                tetrahedron_face_connected_components=int(nc),
                tetrahedron_component_counts=sorted(np.bincount(label).tolist(), reverse=True),
                zero_volume_tetrahedra=zero_count, negative_input_orientation_tetrahedra=signed_negative,
                tetrahedron_faces_with_more_than_two_owners=int(np.count_nonzero(count > 2))), unique[count == 1]


def normal_line_probes(vertices_local_m: np.ndarray, faces: np.ndarray, a: float,
                       b: float, c: float, L: float) -> list[dict]:
    """65 physical infinite-line/triangle intersections, clustered at shared edges.

    No clipping to the analytic interval is allowed.  Entry/exit coordinates
    therefore can lie outside the expected ellipsoid and expose over-thickness.
    """
    points = [(0., 0.)]
    for rho in (.25, .5, .75, .9):
        for j in range(16):
            angle = j*np.pi/8
            points.append((rho*np.cos(angle), rho*np.sin(angle)))
    tri = vertices_local_m[faces]
    origin = tri[:, 0, :2]
    e1 = tri[:, 1, :2]-origin
    e2 = tri[:, 2, :2]-origin
    det = e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0]
    edge_scale = np.linalg.norm(e1, axis=1)*np.linalg.norm(e2, axis=1)
    nonparallel = np.abs(det) > np.maximum(edge_scale*1e-13, L*L*1e-25)
    tol = L*1e-10
    rows = []
    for k, (alpha, beta) in enumerate(points):
        p = np.array([a*alpha, b*beta])
        rr = p-origin
        lam1 = np.zeros(len(faces))
        lam2 = np.zeros(len(faces))
        np.divide(rr[:, 0]*e2[:, 1]-rr[:, 1]*e2[:, 0], det, out=lam1, where=nonparallel)
        np.divide(e1[:, 0]*rr[:, 1]-e1[:, 1]*rr[:, 0], det, out=lam2, where=nonparallel)
        hit = nonparallel & (lam1 >= -1e-9) & (lam2 >= -1e-9) & (lam1+lam2 <= 1+1e-9)
        ss = tri[hit, 0, 2] + lam1[hit]*(tri[hit, 1, 2]-tri[hit, 0, 2]) + lam2[hit]*(tri[hit, 2, 2]-tri[hit, 0, 2])
        ss = np.sort(ss)
        distinct = []
        i = 0
        while i < len(ss):
            j = i+1
            while j < len(ss) and ss[j]-ss[i] <= tol:
                j += 1
            distinct.append(float(np.mean(ss[i:j])))
            i = j
        expected = 2*c*np.sqrt(max(0, 1-alpha*alpha-beta*beta))
        even = len(distinct) % 2 == 0
        actual = sum(distinct[i+1]-distinct[i] for i in range(0, len(distinct)-1, 2)) if even else None
        error = (actual-expected)/expected if actual is not None else None
        entry = distinct[0] if distinct else None
        exit_ = distinct[-1] if distinct else None
        rows.append(dict(probe_id=k, u_over_a=float(alpha), v_over_b=float(beta),
                         rho=float(np.hypot(alpha, beta)), raw_triangle_hits=len(ss),
                         distinct_intersections=len(distinct),
                         actual_interval_count=len(distinct)//2 if even else None,
                         actual_entry_s_m=entry, actual_exit_s_m=exit_,
                         actual_union_thickness_m=actual, expected_thickness_m=float(expected),
                         actual_relative_thickness_error=error,
                         all_intersection_s_m=json.dumps(distinct),
                         intersection_cluster_tolerance_m=tol,
                         pass_one_interval=bool(len(distinct) == 2),
                         pass_thickness_one_percent=bool(error is not None and abs(error) <= .01)))
    return rows


def surface_deviation(vertices_local_m: np.ndarray, faces: np.ndarray, a: float, c: float,
                      output: Path) -> list[dict]:
    """All used vertices, every edge midpoint, centroid and three interior points."""
    barycentric = {"edge_midpoints_with_face_multiplicity": np.array([[.5, .5, 0], [.5, 0, .5], [0, .5, .5]]),
                   "triangle_centroids": np.array([[1/3, 1/3, 1/3]]),
                   "triangle_interior_points": np.array([[.6, .2, .2], [.2, .6, .2], [.2, .2, .6]])}
    results, worst = [], []
    for name in ("vertices", *barycentric):
        count, sum2, maxabs = 0, 0., -1.
        maxpoint, maxsigned = None, None
        if name == "vertices":
            used = np.unique(faces)
            chunks = (vertices_local_m[used[i:i+100_000]] for i in range(0, len(used), 100_000))
        else:
            weights = barycentric[name]
            chunks = (np.einsum("ij,kjl->kil", weights, vertices_local_m[faces[i:i+30_000]]).reshape(-1, 3)
                      for i in range(0, len(faces), 30_000))
        for p in chunks:
            signed = spheroid_signed_distance(p, a, c)
            absolute = np.abs(signed)
            idx = int(np.argmax(absolute))
            if absolute[idx] > maxabs:
                maxabs, maxpoint, maxsigned = float(absolute[idx]), p[idx].tolist(), float(signed[idx])
            count += len(p)
            sum2 += np.dot(signed, signed)
        results.append(dict(sample_category=name, sample_count=count,
                            max_absolute_normal_distance_m=maxabs,
                            rms_normal_distance_m=float(np.sqrt(sum2/count)),
                            max_absolute_distance_over_central_full_thickness=maxabs/(2*c),
                            pass_max_central_thickness_fraction_0p005=bool(maxabs/(2*c) <= .005)))
        worst.append(dict(sample_category=name, local_u_m=maxpoint[0], local_v_m=maxpoint[1],
                          local_n_m=maxpoint[2], signed_nearest_normal_distance_m=maxsigned))
    _csv(output / "geometry_surface_deviation.csv", results)
    _csv(output / "geometry_surface_worst_points_actual.csv", worst)
    return results


def section_segments(vertices_m: np.ndarray, faces: np.ndarray, axis: int, value: float,
                     L: float) -> np.ndarray:
    """Actual plane/triangle intersections; coplanar triangle edges stay visible."""
    q = vertices_m[faces]
    s = q[:, :, axis]-value
    tol = L*1e-11
    chosen = np.flatnonzero((np.min(s, axis=1) <= tol) & (np.max(s, axis=1) >= -tol))
    other = [i for i in range(3) if i != axis]
    segments = []
    for k in chosen:
        verts, d = q[k], s[k]
        hits = []
        if np.max(np.abs(d)) <= tol:
            for i, j in ((0, 1), (1, 2), (2, 0)):
                segments.append(verts[[i, j]][:, other])
            continue
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if abs(d[i]) <= tol:
                hits.append(verts[i])
            if d[i]*d[j] < -tol*tol:
                hits.append(verts[i]+d[i]/(d[i]-d[j])*(verts[j]-verts[i]))
        unique = []
        for p in hits:
            if not any(np.linalg.norm(p-r) <= tol for r in unique):
                unique.append(p)
        if len(unique) == 2:
            segments.append(np.array(unique)[:, other])
    return np.asarray(segments).reshape(-1, 2, 2)


def save_sections(vertices_m: np.ndarray, faces: np.ndarray, spec: dict, output: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    paths = []
    for axis, name in enumerate("xyz"):
        other = [i for i in range(3) if i != axis]
        value = spec["center_m"][axis]
        segments = section_segments(vertices_m, faces, axis, value, spec["L_m"])
        np.savez_compressed(output / f"actual_section_{name}_center.npz", segments_m=segments,
                            fixed_axis=axis, fixed_coordinate_m=value, displayed_axes=other)
        fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
        ax.add_collection(LineCollection(segments*1e6, colors="#182c61", linewidths=.7))
        lo = vertices_m[np.unique(faces)][:, other].min(axis=0)*1e6
        hi = vertices_m[np.unique(faces)][:, other].max(axis=0)*1e6
        margin = max(hi-lo)*.06
        ax.set_xlim(lo[0]-margin, hi[0]+margin)
        ax.set_ylim(lo[1]-margin, hi[1]+margin)
        ax.set_aspect("equal")
        ax.set_xlabel("xyz"[other[0]]+" (um)")
        ax.set_ylabel("xyz"[other[1]]+" (um)")
        scope = "TECHNICAL PILOT, NOT FORMAL CT" if spec["technical_pilot_only"] else "Internal geometry check"
        ax.set_title(f"{scope}: {spec['sample']}, {name}=L/2\nActual exported material-boundary triangle intersections")
        ax.grid(alpha=.2)
        fig.tight_layout()
        file = output / f"actual_section_{name}_center.png"
        fig.savefig(file)
        plt.close(fig)
        paths.append(str(file.resolve()))
    return paths


def verify_geometry(vertices_m: np.ndarray, faces: np.ndarray | None, sample: str,
                    output_dir: str | Path, tetrahedra_inside: np.ndarray | None = None,
                    original_phase_inside: np.ndarray | None = None,
                    metadata: dict | None = None, make_sections: bool = True,
                    parameters: dict | None = None) -> dict:
    """Measure caller-supplied actual geometry.  All arrays use zero-based IDs.

    Exported corner triangles are linear physical chords.  For curved geometric
    elements this is a CHORD CHECK, not a certificate of the high-order surface.
    Supply actual tessellated geometry at controlled accuracy and document it.
    """
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    spec = definition(sample, parameters)
    vertices = _rows(np.asarray(vertices_m, dtype=float), 3, "vertices_m")
    if not np.all(np.isfinite(vertices)):
        raise ValueError("Non-finite vertex")
    tet_stats, boundary_from_tets = None, None
    tets = None
    if tetrahedra_inside is not None:
        tets = _rows(np.asarray(tetrahedra_inside, dtype=np.int64), 4, "tetrahedra_inside")
        if len(tets) == 0 or np.min(tets) < 0 or np.max(tets) >= len(vertices):
            raise ValueError("Empty or invalid inside tetrahedra")
        tet_stats, boundary_from_tets = tetrahedron_checks(vertices, tets)
    if faces is None:
        if boundary_from_tets is None:
            raise ValueError("Supply actual faces or actual inside tetrahedra")
        faces = boundary_from_tets
    faces = _rows(np.asarray(faces, dtype=np.int64), 3, "faces")
    if len(faces) == 0 or np.min(faces) < 0 or np.max(faces) >= len(vertices):
        raise ValueError("Empty or invalid boundary faces")
    topology = face_topology(faces)
    if boundary_from_tets is not None:
        given = np.unique(np.sort(faces, axis=1), axis=0)
        equal = np.array_equal(given, boundary_from_tets)
        topology["surface_equals_inside_tet_union_boundary"] = bool(equal)
    local = (vertices-spec["center_m"]) @ BASIS
    a, b, c = spec["semiaxes_m"]
    probes = normal_line_probes(local, faces, a, b, c, spec["L_m"])
    _csv(output / "geometry_probes_actual.csv", probes)
    deviation = surface_deviation(local, faces, a, c, output)
    # A closed convex boundary is radially oriented relative to its known centre.
    # This reconstruction is disclosed; authoritative tet volume is used if present.
    surface_volume = 0.
    zero_tri = 0
    for start in range(0, len(faces), 100_000):
        q = local[faces[start:start+100_000]]
        cross = np.cross(q[:, 1]-q[:, 0], q[:, 2]-q[:, 0])
        zero_tri += int(np.count_nonzero(np.linalg.norm(cross, axis=1) == 0))
        surface_volume += np.sum(np.abs(np.einsum("ij,ij->i", q[:, 0], cross)))/6
    actual_volume = tet_stats["actual_tetrahedron_volume_m3"] if tet_stats else float(surface_volume)
    used = vertices[np.unique(faces)]
    bbox_min, bbox_max = used.min(axis=0), used.max(axis=0)
    centroid_surface_vertices = used.mean(axis=0)
    local_used = local[np.unique(faces)]
    thickness_errors = [abs(p["actual_relative_thickness_error"]) for p in probes
                        if p["actual_relative_thickness_error"] is not None]
    max_normal = max(d["max_absolute_distance_over_central_full_thickness"] for d in deviation)
    volume_error = actual_volume/spec["ideal_volume_m3"]-1
    all_one = all(p["pass_one_interval"] for p in probes)
    normal_components = tet_stats["tetrahedron_face_connected_components"] if tet_stats else None
    gates = dict(volume_error_le_0p005=abs(volume_error) <= .005,
                 all_65_probes_one_interval=all_one,
                 all_65_thickness_errors_le_0p01=bool(len(thickness_errors) == 65 and max(thickness_errors) <= .01),
                 sampled_normal_error_over_full_thickness_le_0p005=max_normal <= .005,
                 closed_connected_surface=bool(topology["closed_two_manifold"] and topology["surface_edge_connected_components"] == 1 and zero_tri == 0),
                 union_tetrahedra_face_connected=None if tet_stats is None else normal_components == 1,
                 wholly_within_central_half_z=bool(bbox_min[2] > .25*spec["L_m"] and bbox_max[2] < .75*spec["L_m"]),
                 wholly_within_protected_core=(bool(bbox_min[2] > .25*spec["L_m"] and bbox_max[2] < .75*spec["L_m"])
                                              if spec["protection_regions_present"] else None))
    if boundary_from_tets is not None:
        gates["boundary_matches_actual_inside_tetrahedra"] = topology["surface_equals_inside_tet_union_boundary"]
    summary = dict(schema="r7b_actual_boundary_geometry_v1", sample=sample,
                   status="ACTUAL_INPUT_MESH_MEASURED_NOT_A_COMSOL_SOLVE",
                   geometry_scope=spec["geometry_scope"], technical_pilot_only=spec["technical_pilot_only"],
                   source_metadata=metadata or {}, definition=spec,
                   input_geometry_representation="linear exported triangle/tetrahedron corners",
                   volume_measurement_method="sum absolute inside tetrahedron volumes" if tet_stats else "closed convex surface, faces radially oriented about prescribed centre",
                   actual_volume_m3=actual_volume, ideal_volume_m3=spec["ideal_volume_m3"],
                   actual_volume_relative_error=volume_error,
                   actual_radially_oriented_surface_volume_m3=float(surface_volume),
                   actual_bbox_min_m=bbox_min, actual_bbox_max_m=bbox_max,
                   actual_bbox_min_error_m=bbox_min-spec["ideal_bbox_min_m"],
                   actual_bbox_max_error_m=bbox_max-spec["ideal_bbox_max_m"],
                   actual_prescribed_axis_projection_min_m=local_used.min(axis=0),
                   actual_prescribed_axis_projection_max_m=local_used.max(axis=0),
                   actual_prescribed_axis_extent_m=np.ptp(local_used, axis=0),
                   surface_vertex_centroid_m=centroid_surface_vertices,
                   surface_vertex_centroid_is_not_volume_centroid=True,
                   actual_central_thickness_m=probes[0]["actual_union_thickness_m"],
                   actual_central_thickness_relative_error=probes[0]["actual_relative_thickness_error"],
                   maximum_probe_absolute_relative_thickness_error=max(thickness_errors) if thickness_errors else None,
                   maximum_sampled_normal_error_over_central_thickness=max_normal,
                   zero_area_boundary_triangles=zero_tri, topology=topology,
                   inside_tetrahedron_checks=tet_stats, gates=gates,
                   all_evaluated_geometric_gates_pass=all(v for v in gates.values() if v is not None),
                   tetrahedron_connectivity_independently_checked=tet_stats is not None,
                   limitations=["Sampling the boundary is not a rigorous global Hausdorff bound.",
                                "Closed face topology and probes do not alone exclude arbitrary triangle self-intersections.",
                                "This helper does not certify preservation of the external original CT phase field.",
                                "Geometry checks do not establish electrical convergence.",
                                "With curved elements, corner-chord measurements need separate interpretation."])
    if original_phase_inside is not None:
        if tets is None:
            raise ValueError("original_phase_inside requires tetrahedra_inside")
        labels = np.asarray(original_phase_inside).ravel()
        if len(labels) != len(tets) or not np.all(np.isin(labels, [0, 1])):
            raise ValueError("original_phase_inside must contain one 0/1 original label per inside tet")
        vols = {"original_pore_m3": 0., "original_skeleton_m3": 0.}
        for start in range(0, len(tets), 100_000):
            q = vertices[tets[start:start+100_000]]
            vv = np.abs(np.einsum("ij,ij->i", q[:, 1]-q[:, 0], np.cross(q[:, 2]-q[:, 0], q[:, 3]-q[:, 0])))/6
            lab = labels[start:start+100_000]
            vols["original_pore_m3"] += float(vv[lab == 0].sum())
            vols["original_skeleton_m3"] += float(vv[lab == 1].sum())
        summary["actual_ellipsoid_original_phase_overlap"] = vols
    summary["internal_section_images"] = save_sections(vertices, faces, spec, output) if make_sections else []
    _dump(output / "geometry_verification.json", summary)
    row = dict(sample=sample, geometry_scope=spec["geometry_scope"],
               actual_sample_L_m=spec["L_m"], actual_voxel_size_m=spec["voxel_size_m"],
               actual_volume_m3=actual_volume,
               ideal_volume_m3=spec["ideal_volume_m3"], relative_volume_error=volume_error,
               actual_central_thickness_m=summary["actual_central_thickness_m"],
               actual_central_thickness_relative_error=summary["actual_central_thickness_relative_error"],
               maximum_probe_absolute_relative_thickness_error=summary["maximum_probe_absolute_relative_thickness_error"],
               maximum_sampled_normal_error_over_central_thickness=max_normal,
               surface_components=topology["surface_edge_connected_components"],
               tetrahedron_face_components=normal_components,
               all_evaluated_geometric_gates_pass=summary["all_evaluated_geometric_gates_pass"],
               original_phase_preservation="not established by this helper")
    _csv(output / "geometry_verification.csv", [row])
    return summary


def load_input(path: Path, index_base: int | None = None) -> dict:
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as z:
            values = {k: z[k] for k in z.files}
    elif path.suffix.lower() == ".mat":
        from scipy.io import loadmat
        try:
            values = loadmat(path)
        except NotImplementedError:
            raise ValueError("Use MATLAB save -v7 for the small boundary export; v7.3 requires h5py, not this environment")
    else:
        raise ValueError("Input must be .npz or MATLAB -v7 .mat")
    base = int(np.asarray(values.get("index_base", 0 if index_base is None else index_base)).squeeze())
    if index_base is not None:
        base = index_base
    if base not in (0, 1):
        raise ValueError("index_base must be 0 or 1")
    v = values.get("vertices_m")
    if v is None:
        raise ValueError("vertices_m key required; automatic unit guessing is forbidden")
    kwargs = dict(vertices_m=v, faces=None)
    if "faces" in values:
        kwargs["faces"] = np.asarray(values["faces"], dtype=np.int64)-base
    if "tetrahedra_inside" in values:
        kwargs["tetrahedra_inside"] = np.asarray(values["tetrahedra_inside"], dtype=np.int64)-base
    if "original_phase_inside" in values:
        kwargs["original_phase_inside"] = values["original_phase_inside"]
    return kwargs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--sample", choices=SPACING, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--index-base", type=int, choices=(0, 1))
    parser.add_argument("--no-sections", action="store_true")
    parser.add_argument("--parameters", type=Path,
                        help="Actual source parameters.json (L_m, spacing_m, optional original shape/axes); required for technical32")
    args = parser.parse_args()
    parameters = json.loads(args.parameters.read_text(encoding="utf-8-sig")) if args.parameters else None
    result = verify_geometry(**load_input(args.input, args.index_base), sample=args.sample,
                             output_dir=args.output, make_sections=not args.no_sections,
                             parameters=parameters,
                             metadata={"actual_mesh_export": str(args.input.resolve()),
                                       "actual_parameters_file": str(args.parameters.resolve()) if args.parameters else None})
    print(json.dumps({k: result[k] for k in ("sample", "geometry_scope", "actual_volume_relative_error",
                     "actual_central_thickness_relative_error", "maximum_probe_absolute_relative_thickness_error",
                     "maximum_sampled_normal_error_over_central_thickness", "all_evaluated_geometric_gates_pass")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

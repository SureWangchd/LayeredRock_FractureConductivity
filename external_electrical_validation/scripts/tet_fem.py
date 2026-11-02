"""Linear tetrahedron FEM for div(sigma grad(u))=0, natural insulation elsewhere.

Coordinates are in metres and conductivities in S/m when SI conductance is needed.
``energy`` is u.T K u = integral grad(u).sigma.grad(u) dV (electrical power),
not half of that value. Dirichlet reactions are K u at fixed nodes; the sum on
the high-potential electrode is positive current entering the domain.

No elements are removed or repaired. Degenerate elements and invalid materials
raise errors. Arbitrary vertex orientation is supported through |det(J)|.
"""

from __future__ import annotations

import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg, spsolve, LinearOperator


def _geometry(nodes, cells):
    nodes = np.asarray(nodes, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or not np.isfinite(nodes).all():
        raise ValueError("nodes must be a finite N x 3 array")
    if cells.ndim != 2 or cells.shape[1] != 4 or cells.size == 0:
        raise ValueError("cells must be a nonempty M x 4 array")
    if cells.min() < 0 or cells.max() >= len(nodes):
        raise ValueError("cell indices are outside nodes")
    return nodes, cells


def _material(conductivity, ncells):
    a = np.asarray(conductivity, dtype=np.float64)
    if a.ndim == 0:
        if not np.isfinite(a) or a <= 0:
            raise ValueError("conductivity must be finite and strictly positive")
        return "scalar", float(a)
    if a.shape == (ncells,):
        if not np.isfinite(a).all() or (a <= 0).any():
            raise ValueError("all cell conductivities must be finite and positive")
        return "cell_scalar", a
    if a.shape not in ((3, 3), (ncells, 3, 3)):
        raise ValueError("conductivity must be scalar, M, 3x3, or Mx3x3")
    if not np.isfinite(a).all():
        raise ValueError("non-finite conductivity tensor")
    if not np.allclose(a, np.swapaxes(a, -1, -2), rtol=1e-12, atol=0):
        raise ValueError("conductivity tensors must be symmetric")
    try:
        np.linalg.cholesky(a)
    except np.linalg.LinAlgError as exc:
        raise ValueError("conductivity tensors must be positive definite") from exc
    return ("tensor" if a.ndim == 2 else "cell_tensor"), a


def element_gradients(nodes, cells):
    """Return gradients (M,4,3), positive volumes (M,), signed det(J) (M,).

    J[:,j] = x[j+1]-x[0].  The reference gradient of N[j+1] is e_j,
    so its physical *row* gradient is row j of inv(J), not inv(J).T.
    """
    x = np.asarray(nodes)[np.asarray(cells)]
    jac = np.swapaxes(x[:, 1:, :] - x[:, :1, :], 1, 2)
    det = np.linalg.det(jac)
    if not np.isfinite(det).all() or (det == 0).any():
        bad = np.flatnonzero(~np.isfinite(det) | (det == 0))
        raise ValueError(f"zero-volume/nonfinite tetrahedra: {bad[:20].tolist()}")
    inv = np.linalg.inv(jac)
    grad = np.empty((len(x), 4, 3), dtype=np.float64)
    grad[:, 1:, :] = inv
    grad[:, 0, :] = -inv.sum(axis=1)
    return grad, np.abs(det) / 6.0, det


def assemble(nodes, cells, conductivity, *, chunk_size=100_000,
             return_info=False):
    """Assemble symmetric CSR stiffness matrix; all supplied tetrahedra are used."""
    nodes, cells = _geometry(nodes, cells)
    kind, sigma = _material(conductivity, len(cells))
    n = len(nodes)
    # The full COO storage is 16 entries per tetrahedron; use 32-bit indices when
    # possible. Local geometry is evaluated in chunks to bound temporary memory.
    idx_type = np.int32 if n <= np.iinfo(np.int32).max else np.int64
    rows = np.empty(16 * len(cells), dtype=idx_type)
    cols = np.empty_like(rows)
    data = np.empty(16 * len(cells), dtype=np.float64)
    min_volume, max_volume, volume, reversed_count = np.inf, 0.0, 0.0, 0
    for start in range(0, len(cells), int(chunk_size)):
        stop = min(start + int(chunk_size), len(cells))
        c = cells[start:stop]
        grad, vol, det = element_gradients(nodes, c)
        if kind == "scalar":
            local = np.einsum("mik,mjk->mij", grad, grad) * sigma
        elif kind == "cell_scalar":
            local = np.einsum("mik,mjk->mij", grad, grad) * sigma[start:stop, None, None]
        elif kind == "tensor":
            local = np.einsum("mik,kl,mjl->mij", grad, sigma, grad, optimize=True)
        else:
            local = np.einsum("mik,mkl,mjl->mij", grad, sigma[start:stop], grad, optimize=True)
        local *= vol[:, None, None]
        sl = slice(16 * start, 16 * stop)
        rows[sl] = np.repeat(c, 4, axis=1).reshape(-1)
        cols[sl] = np.tile(c, (1, 4)).reshape(-1)
        data[sl] = local.reshape(-1)
        min_volume = min(min_volume, float(vol.min()))
        max_volume = max(max_volume, float(vol.max()))
        volume += float(vol.sum())
        reversed_count += int((det < 0).sum())
    matrix = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    info = dict(nodes=n, elements=len(cells), volume=volume,
                min_element_volume=min_volume, max_element_volume=max_volume,
                negative_orientation_count=reversed_count, nnz=matrix.nnz,
                conductivity_kind=kind)
    return (matrix, info) if return_info else matrix


def _constraints(fixed_nodes, fixed_values, n):
    ids = np.asarray(fixed_nodes, dtype=np.int64).reshape(-1)
    vals = np.asarray(fixed_values, dtype=np.float64)
    if vals.ndim == 0:
        vals = np.full(len(ids), float(vals))
    vals = vals.reshape(-1)
    if len(ids) == 0 or len(ids) != len(vals):
        raise ValueError("nonempty Dirichlet nodes/values with matching lengths required")
    if ids.min() < 0 or ids.max() >= n or not np.isfinite(vals).all():
        raise ValueError("invalid Dirichlet index/value")
    order = np.argsort(ids)
    ids, vals = ids[order], vals[order]
    duplicate = ids[1:] == ids[:-1]
    if np.any(duplicate & (vals[1:] != vals[:-1])):
        raise ValueError("conflicting Dirichlet values at a node")
    unique = np.r_[True, ~duplicate]
    return ids[unique], vals[unique]


def solve(nodes, cells, conductivity, fixed_nodes, fixed_values, *,
          method="auto", rtol=1e-10, atol=0.0, maxiter=5000,
          direct_threshold=80_000, use_pyamg=True, chunk_size=100_000,
          return_matrix=False):
    """Solve the electrical Dirichlet problem, returning potentials and diagnostics.

    Optional method: 'auto', 'direct', or 'cg'. CG uses pyamg when available,
    otherwise diagonal scaling. Nonconvergence raises RuntimeError; it cannot
    silently be reported as a successful prediction. ``relative_residual`` is
    ||K_ff u_f-b_f||_2 / ||b_f||_2 (denominator 1 for zero b_f).
    """
    t0 = time.perf_counter()
    matrix, mesh_info = assemble(nodes, cells, conductivity,
                                 chunk_size=chunk_size, return_info=True)
    n = matrix.shape[0]
    ids, vals = _constraints(fixed_nodes, fixed_values, n)
    free_mask = np.ones(n, dtype=bool)
    free_mask[ids] = False
    free = np.flatnonzero(free_mask)
    u = np.zeros(n)
    u[ids] = vals
    # Setting boundary values in a full vector avoids materializing K_fc.
    rhs = -np.asarray(matrix @ u)[free]
    kff = matrix[free][:, free].tocsr()
    if method not in ("auto", "direct", "cg"):
        raise ValueError("method must be auto, direct, or cg")
    selected = ("direct" if len(free) <= direct_threshold else "cg") if method == "auto" else method
    iterations, preconditioner = 0, "none"
    if len(free):
        if selected == "direct":
            x = spsolve(kff, rhs)
        else:
            diag = kff.diagonal()
            if np.any(diag <= 0):
                raise ValueError("nonpositive diagonal: unused/disconnected mesh nodes")
            prec = None
            if use_pyamg:
                try:
                    import pyamg
                except ImportError:
                    pass
                else:
                    hierarchy = pyamg.smoothed_aggregation_solver(kff, symmetry="symmetric")
                    prec = hierarchy.aspreconditioner(cycle="V")
                    preconditioner = "pyamg_smoothed_aggregation"
            if prec is None:
                prec = LinearOperator(kff.shape, matvec=lambda y: y / diag)
                preconditioner = "jacobi"
            def count(_):
                nonlocal iterations
                iterations += 1
            try:
                x, status = cg(kff, rhs, M=prec, rtol=rtol, atol=atol,
                               maxiter=maxiter, callback=count)
            except TypeError:  # scipy before the rtol keyword transition
                x, status = cg(kff, rhs, M=prec, tol=rtol, atol=atol,
                               maxiter=maxiter, callback=count)
            if status != 0:
                raise RuntimeError(f"CG did not converge: status={status}, iterations={iterations}")
        if not np.isfinite(x).all():
            raise RuntimeError("nonfinite FEM solution; check mesh connectivity and constraints")
        u[free] = x
        residual = kff @ x - rhs
        relative_residual = float(np.linalg.norm(residual) / (np.linalg.norm(rhs) or 1.0))
        absolute_residual = float(np.linalg.norm(residual))
    else:
        relative_residual = absolute_residual = 0.0
    reaction = np.asarray(matrix @ u)
    energy = float(u @ reaction)
    boundary_power = float(vals @ reaction[ids])
    result = dict(u=u, fixed_nodes=ids, fixed_values=vals,
                  dirichlet_reaction=reaction[ids], reaction=reaction,
                  energy=energy, half_energy=0.5 * energy,
                  boundary_power=boundary_power,
                  power_balance_absolute=abs(energy - boundary_power),
                  relative_residual=relative_residual,
                  absolute_residual=absolute_residual,
                  current_balance_absolute=abs(float(reaction[ids].sum())),
                  solver_info=dict(method=selected, preconditioner=preconditioner,
                                   iterations=iterations, rtol=rtol, atol=atol,
                                   seconds=time.perf_counter()-t0),
                  mesh_info=mesh_info)
    if return_matrix:
        result["matrix"] = matrix
    return result


def _transmission_check():
    """Two fitted materials: exact series and parallel conduction, no averaging.

    The rectangular grid has an existing plane at each interface. Element phase
    assignment selects whole tetrahedra that already lie on one side of that
    plane; no tetrahedron crosses a material interface.
    """
    import itertools
    lengths = np.array([0.021, 0.013, 0.008])
    divisions = [6, 4, 4]
    grid = np.meshgrid(*(np.linspace(0, l, k+1) for l, k in zip(lengths, divisions)), indexing="ij")
    nodes = np.stack(grid, axis=-1).reshape(-1, 3)
    shape = tuple(k+1 for k in divisions)
    cells = []
    for base in itertools.product(*(range(k) for k in divisions)):
        for perm in itertools.permutations(range(3)):
            vertex = np.array(base)
            tet = [np.ravel_multi_index(vertex, shape)]
            for axis in perm:
                vertex = vertex.copy()
                vertex[axis] += 1
                tet.append(np.ravel_multi_index(vertex, shape))
            cells.append(tet)
    cells = np.asarray(cells)
    low = np.flatnonzero(nodes[:, 0] == 0)
    high = np.flatnonzero(nodes[:, 0] == lengths[0])
    fixed = np.r_[low, high]
    values = np.r_[np.zeros(len(low)), np.ones(len(high))]
    sigma0, sigma1 = 0.28, 4.69
    gradient, volumes, _ = element_gradients(nodes, cells)
    records = []
    for name, interface_axis in [("series_interface_transmission", 0), ("parallel_interface_transmission", 2)]:
        position = lengths[interface_axis] / 2
        vertices = nodes[cells, interface_axis]
        # Strict open-side check: the interface itself is shared by both phases.
        crosses = (vertices.min(axis=1) < position-1e-14) & (vertices.max(axis=1) > position+1e-14)
        assert not crosses.any(), "A tetrahedron crosses the material interface"
        material_one = vertices.mean(axis=1) > position
        sigma = np.where(material_one, sigma1, sigma0)
        out = solve(nodes, cells, sigma, fixed, values, method="direct")
        if interface_axis == 0:
            equivalent_sigma = 1/(0.5/sigma0+0.5/sigma1)
            resistance_length = lengths[0]*(0.5/sigma0+0.5/sigma1)
            expected_u = (np.minimum(nodes[:,0],position)/sigma0 + np.maximum(nodes[:,0]-position,0)/sigma1)/resistance_length
            expected_flux_x = np.full(len(cells), -equivalent_sigma/lengths[0])
        else:
            equivalent_sigma = 0.5*(sigma0+sigma1)
            expected_u = nodes[:,0]/lengths[0]
            expected_flux_x = -sigma/lengths[0]
        conductance = float(out["dirichlet_reaction"][out["fixed_values"] == 1].sum())
        exact_conductance = equivalent_sigma*lengths[1]*lengths[2]/lengths[0]
        grad_u = np.einsum("mi,mij->mj", out["u"][cells], gradient)
        current_density = -sigma[:,None]*grad_u
        flux_error = float(abs(current_density[:,0]-expected_flux_x).max())
        transverse_error = float(abs(current_density[:,1:]).max())
        relative_g_error = abs(conductance/exact_conductance-1)
        potential_error = float(abs(out["u"]-expected_u).max())
        assert relative_g_error < 1e-10, (name, relative_g_error)
        assert potential_error < 1e-11, (name, potential_error)
        assert flux_error < 1e-8, (name, flux_error)
        assert transverse_error < 1e-8, (name, transverse_error)
        records.append(dict(name=name,conductivities_S_per_m=[sigma0,sigma1],
                            phase_volume_fractions=[float(volumes[~material_one].sum()/volumes.sum()),float(volumes[material_one].sum()/volumes.sum())],
                            analytic_effective_conductivity_S_per_m=float(equivalent_sigma),
                            analytic_conductance_S=float(exact_conductance),numerical_conductance_S=conductance,
                            relative_conductance_error=float(relative_g_error),maximum_potential_error=potential_error,
                            maximum_flux_error_A_per_m2=flux_error,maximum_transverse_current_A_per_m2=transverse_error,
                            interface_crossing_elements=int(crosses.sum()),relative_residual=out["relative_residual"],
                            energy_W=out["energy"],boundary_power_W=out["boundary_power"],
                            power_balance_absolute_W=out["power_balance_absolute"],mesh=out["mesh_info"]))
    return records


def _self_check():
    """Analytic box conductance and affine skew-tetrahedron gradient checks."""
    import itertools
    lengths = np.array([0.021, 0.013, 0.008])
    divisions = [4, 3, 2]
    grid = np.meshgrid(*(np.linspace(0, l, k+1) for l, k in zip(lengths, divisions)), indexing="ij")
    nodes = np.stack(grid, axis=-1).reshape(-1, 3)
    shape = tuple(k+1 for k in divisions)
    cells = []
    for base in itertools.product(*(range(k) for k in divisions)):
        for perm in itertools.permutations(range(3)):
            vertex = np.array(base)
            tet = [np.ravel_multi_index(vertex, shape)]
            for axis in perm:
                vertex = vertex.copy()
                vertex[axis] += 1
                tet.append(np.ravel_multi_index(vertex, shape))
            cells.append(tet)
    cells = np.asarray(cells)
    low = np.flatnonzero(nodes[:, 0] == 0)
    high = np.flatnonzero(nodes[:, 0] == lengths[0])
    fixed = np.r_[low, high]
    values = np.r_[np.zeros(len(low)), np.ones(len(high))]
    expected_u = nodes[:, 0] / lengths[0]
    records = []
    for name, conductivity, sigmax, method in [
        ("homogeneous_scalar_direct", 0.7, 0.7, "direct"),
        ("homogeneous_diagonal_tensor_direct", np.diag([0.7, 0.31, 1.2]), 0.7, "direct"),
        ("homogeneous_scalar_cg", np.full(len(cells), 0.7), 0.7, "cg"),
        ("cell_tensor_direct", np.broadcast_to(np.diag([0.7,0.31,1.2]),(len(cells),3,3)), 0.7, "direct")]:
        out = solve(nodes, cells, conductivity, fixed, values, method=method, use_pyamg=False)
        conductance = float(out["dirichlet_reaction"][out["fixed_values"] == 1].sum())
        analytic = sigmax * lengths[1] * lengths[2] / lengths[0]
        rec = dict(name=name, numerical_conductance_S=conductance,
                   analytic_conductance_S=float(analytic),
                   relative_conductance_error=abs(conductance/analytic-1),
                   maximum_potential_error=float(abs(out["u"]-expected_u).max()),
                   relative_residual=out["relative_residual"],
                   power_error=abs(out["energy"]-analytic),
                   mesh=out["mesh_info"])
        assert rec["relative_conductance_error"] < 1e-10, rec
        assert rec["maximum_potential_error"] < 1e-10, rec
        records.append(rec)
    # Nonorthogonal affine map makes a mistaken inverse transpose detectable.
    skew = np.array([[0.2,-0.4,0.1],[1.2,0.1,0.3],[0.4,0.7,-0.2],[-0.1,-0.2,1.3]])
    tc = np.array([[0,1,2,3],[0,2,1,3]])
    gradient, volume, det = element_gradients(skew, tc)
    exact_gradient = np.array([1.7,-2.3,0.6])
    u = skew @ exact_gradient + 0.4
    reconstructed = np.einsum("mi,mij->mj", u[tc], gradient)
    error = float(abs(reconstructed-exact_gradient).max())
    assert error < 1e-13, error
    tensor = np.array([[2.,0.3,0.1],[0.3,1.2,-0.2],[0.1,-0.2,0.9]])
    k = assemble(skew, tc[:1], tensor)
    expected_energy = float(volume[0] * (exact_gradient @ tensor @ exact_gradient))
    energy_error = abs(float(u @ k @ u)-expected_energy)
    assert energy_error < 1e-13, energy_error
    records.append(dict(name="skew_tetra_gradient_and_tensor_energy",
                        maximum_gradient_error=error, energy_error=energy_error,
                        signed_determinants=det.tolist(), volumes=volume.tolist()))
    records.extend(_transmission_check())
    return dict(passed=True, checks=records,
                energy_definition="integral grad(u).sigma.grad(u) dV, no 1/2",
                notes="Full-face electrodes on uniform boxes; affine patch test on skew tetrahedra including reversed vertex order; fitted two-phase interface series/parallel analytical transmission tests.")


if __name__ == "__main__":
    import json
    import sys
    from pathlib import Path
    result = _self_check()
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))

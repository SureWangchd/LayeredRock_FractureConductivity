"""Layered-exterior response helpers for the proposed fifth research round.

All dimensional inputs are SI. This module DOES NOT run COMSOL or reproduce
any rock's full electric response. It evaluates analytical exterior coefficients
and illustrates the algebra used to couple them to a supplied core operator.
Layers are ordered from the CORE interface towards the fixed-potential ELECTRODE.

For mode n, psi_n(x)=cos(n*pi*x/W), q_n=Y_n*v_n (n>0).
For n=0, q_0=Y_0*(v_0-V_electrode).
q is outward current density from the core into the exterior [A/m^2].
Y has units S/m^2. Core matrices below already include thickness [m].
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import math
import numpy as np


@dataclass(frozen=True)
class Layer:
    thickness_m: float
    conductivity_S_m: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.thickness_m) or self.thickness_m <= 0:
            raise ValueError("Layer thickness must be positive and finite [m].")
        if not math.isfinite(self.conductivity_S_m) or self.conductivity_S_m <= 0:
            raise ValueError("Conductivity must be positive and finite [S/m].")


def _layers(layers: Sequence[Layer]) -> tuple[Layer, ...]:
    result = tuple(layers)
    if not result or not all(isinstance(x, Layer) for x in result):
        raise ValueError("Provide a nonempty sequence of Layer objects.")
    return result


def input_admittance(layers: Sequence[Layer], k_per_m: float) -> float:
    """Exact scalar DtN coefficient for flat isotropic layers.

    The far electrode has ZERO modal perturbation. A constant electrode voltage
    must be handled by the affine n=0 term, NOT by a nonzero-mode source.
    Uses a stable tanh recursion, not products of cosh transfer matrices.
    """
    layers = _layers(layers)
    if not math.isfinite(k_per_m) or k_per_m < 0:
        raise ValueError("k must be nonnegative and finite [1/m].")
    if k_per_m == 0:
        return 1.0 / sum(a.thickness_m / a.conductivity_S_m for a in layers)
    y = None  # Dirichlet load at the far electrode: admittance -> infinity.
    for layer in reversed(layers):
        a = layer.conductivity_S_m * k_per_m
        th = math.tanh(k_per_m * layer.thickness_m)
        if y is None:
            y = a / th
        else:
            # Algebraically a*(y+a*th)/(a+y*th), scaled for stability.
            scale = max(a, y)
            aa, yy = a / scale, y / scale
            y = a * (yy + aa * th) / (aa + yy * th)
    return float(y)


def cosine_norms(width_m: float, nmax: int) -> np.ndarray:
    if not math.isfinite(width_m) or width_m <= 0:
        raise ValueError("width_m must be positive and finite.")
    if not isinstance(nmax, int) or isinstance(nmax, bool) or nmax < 0:
        raise ValueError("nmax must be a nonnegative integer.")
    n = np.full(nmax + 1, width_m / 2.0)
    n[0] = width_m
    return n


def modal_admittances(layers: Sequence[Layer], width_m: float,
                       nmax: int) -> np.ndarray:
    cosine_norms(width_m, nmax)  # Validation.
    return np.array([input_admittance(layers, n * math.pi / width_m)
                     for n in range(nmax + 1)])


def equivalent_tensor(layers: Sequence[Layer]) -> tuple[float, float, float]:
    """Return H [m], sigma_x and sigma_z [S/m], arithmetic/harmonic means."""
    layers = _layers(layers)
    h = sum(a.thickness_m for a in layers)
    sx = sum(a.thickness_m * a.conductivity_S_m for a in layers) / h
    sz = h / sum(a.thickness_m / a.conductivity_S_m for a in layers)
    return h, sx, sz


def equivalent_modal_admittances(layers: Sequence[Layer], width_m: float,
                                  nmax: int) -> np.ndarray:
    """Common homogeneous ANISOTROPIC exterior, preserving both DC directions."""
    cosine_norms(width_m, nmax)
    h, sx, sz = equivalent_tensor(layers)
    y = np.empty(nmax + 1)
    y[0] = sz / h
    for n in range(1, nmax + 1):
        k = n * math.pi / width_m
        alpha = k * math.sqrt(sx / sz)
        y[n] = sz * alpha / math.tanh(alpha * h)
    return y


def hybrid_admittances(exact: np.ndarray, reference: np.ndarray,
                        retained_order: int) -> np.ndarray:
    """Use exact n=0..K and a COMMON reference for the rest of this array.

    IMPORTANT: modes ABOVE THE ARRAY END must also retain a reference operator
    in the production FEM model. This array is NOT permission to set the tail
    to zero or to make it insulating. See 02_math_and_implementation_CN.md.
    """
    exact, reference = np.asarray(exact, float), np.asarray(reference, float)
    if exact.ndim != 1 or exact.shape != reference.shape or exact.size < 1:
        raise ValueError("Admittance arrays must have the same 1D shape.")
    if np.any(~np.isfinite(exact)) or np.any(~np.isfinite(reference)):
        raise ValueError("Admittances must be finite.")
    if np.any(exact <= 0) or np.any(reference <= 0):
        raise ValueError("Both full admittance maps must be positive.")
    if not isinstance(retained_order, int) or not 0 <= retained_order < exact.size:
        raise ValueError("retained_order must be in [0, len(exact)-1].")
    out = reference.copy()
    out[:retained_order + 1] = exact[:retained_order + 1]
    return out


def modal_update_weights(delta_y: np.ndarray, width_m: float,
                          thickness_m: float) -> np.ndarray:
    """Return weights w_n=t*DeltaY_n/N_n for the FEM correction C diag(w) C.T.

    C[i,n]=integral_Gamma N_i(x)*cos(n*pi*x/W) dx, in METRES.
    DeltaY may be negative; clipping it changes the physics.
    This only builds a correction. A positive reference exterior must remain.
    """
    dy = np.asarray(delta_y, float)
    if dy.ndim != 1 or dy.size == 0 or not np.all(np.isfinite(dy)):
        raise ValueError("delta_y must be a nonempty finite 1D array.")
    if not math.isfinite(thickness_m) or thickness_m <= 0:
        raise ValueError("thickness_m must be positive and finite.")
    return thickness_m * dy / cosine_norms(width_m, dy.size - 1)


def apply_modal_update(C: np.ndarray, weights: np.ndarray,
                        potential: np.ndarray) -> np.ndarray:
    """Apply C diag(weights) C.T without allocating a dense global matrix."""
    C, w, u = np.asarray(C, float), np.asarray(weights, float), np.asarray(potential, float)
    if C.ndim != 2 or w.shape != (C.shape[1],) or u.shape != (C.shape[0],):
        raise ValueError("Expected C:(ndof,nmodes), weights:(nmodes,), u:(ndof,).")
    if not (np.all(np.isfinite(C)) and np.all(np.isfinite(w)) and np.all(np.isfinite(u))):
        raise ValueError("All matrix/vector entries must be finite.")
    return C @ (w * (C.T @ u))


def band_layers(distance_m: float, host_S_m: float, band_S_m: float,
                band_thickness_m: float = 0.002,
                exterior_height_m: float = 0.030) -> tuple[Layer, ...]:
    """An isolated band in a homogeneous host, ordered core -> electrode."""
    if not (0 < distance_m < exterior_height_m - band_thickness_m):
        raise ValueError("Band must lie strictly inside the exterior.")
    return (Layer(distance_m, host_S_m), Layer(band_thickness_m, band_S_m),
            Layer(exterior_height_m - distance_m - band_thickness_m, host_S_m))


def original_transfer_check(layers: Sequence[Layer], k_per_m: float) -> float:
    """Independent transfer-matrix check, suitable ONLY for small/moderate kH.

    State is [V, p] with p=sigma*dV/ds, s pointing core -> electrode.
    M maps the core state to the electrode; V_e=0 implies Y=M00/M01.
    Large-k production runs must use input_admittance instead.
    """
    layers = _layers(layers)
    if k_per_m < 0 or not math.isfinite(k_per_m):
        raise ValueError("Invalid k.")
    if k_per_m * sum(a.thickness_m for a in layers) > 15:
        raise ValueError("Transfer check restricted to kH<=15.")
    M = np.eye(2)
    for a in layers:
        if k_per_m == 0:
            T = np.array([[1.0, a.thickness_m / a.conductivity_S_m], [0., 1.]])
        else:
            kh = k_per_m * a.thickness_m
            c, s = math.cosh(kh), math.sinh(kh)
            ak = a.conductivity_S_m * k_per_m
            T = np.array([[c, s / ak], [ak * s, c]])
        M = T @ M
    return float(M[0, 0] / M[0, 1])


def solve_supplied_modal_core(core_S: np.ndarray, lower_Y: np.ndarray,
                              upper_Y: np.ndarray, width_m: float,
                              thickness_m: float = 1.0, lower_V: float = 0.0,
                              upper_V: float = 1.0) -> dict:
    """Algebraic demonstration for an ALREADY CONDENSED modal core operator.

    core_S is a symmetric PSD Dirichlet-energy matrix [S] for modal trace
    coefficients ordered [lower n=0..N, upper n=0..N], INCLUDING thickness.
    It must come from an independent core calculation, never a target field.
    A finite modal trace space suppresses higher traces; this helper does not
    implement the recommended full-reference/high-mode-completion production
    method. Use it for algebra checks or an explicitly labeled Galerkin model.
    """
    lo, hi, S = np.asarray(lower_Y, float), np.asarray(upper_Y, float), np.asarray(core_S, float)
    if lo.ndim != 1 or hi.shape != lo.shape or lo.size == 0:
        raise ValueError("lower_Y and upper_Y must be equal-length vectors.")
    m = lo.size
    if S.shape != (2*m, 2*m):
        raise ValueError("core_S shape must be (2*nmodes,2*nmodes).")
    if np.any(lo <= 0) or np.any(hi <= 0) or not np.all(np.isfinite(S)):
        raise ValueError("Positive admittances and finite core matrix required.")
    if not np.allclose(S, S.T, atol=1e-12, rtol=1e-10):
        raise ValueError("core_S must be symmetric.")
    scale = max(float(np.linalg.norm(S, 2)), 1.0)
    if float(np.linalg.eigvalsh(S).min()) < -1e-11 * scale:
        raise ValueError("core_S must be positive semidefinite.")
    if thickness_m <= 0 or upper_V == lower_V:
        raise ValueError("Positive thickness and nonzero applied voltage required.")
    norms = cosine_norms(width_m, m-1)
    weights = thickness_m * np.concatenate((norms*lo, norms*hi))
    electrode = np.zeros(2*m)
    electrode[0], electrode[m] = lower_V, upper_V
    A = S + np.diag(weights)
    b = weights * electrode
    v = np.linalg.solve(A, b)
    residual = A @ v - b
    dv = upper_V - lower_V
    pcore = float(v @ S @ v)
    pext = float(np.dot(weights, (v-electrode)**2))
    iupper = thickness_m*width_m*hi[0]*(upper_V-v[m])
    ilower = thickness_m*width_m*lo[0]*(v[0]-lower_V)
    return {"trace_coefficients_V": v, "G_S": (pcore+pext)/dv**2,
            "core_power_W": pcore, "exterior_power_W": pext,
            "upper_inward_current_A": iupper, "lower_outward_current_A": ilower,
            "linear_residual_norm_A": float(np.linalg.norm(residual))}


if __name__ == "__main__":
    import json
    W, nmax = 0.120, 8
    out = {"scope": "Analytical exterior only; NOT new COMSOL/full-rock results.",
           "mode_definition": "cos(n*pi*x/W), W=0.120 m; Y in S/m^2", "sides": {}}
    for side, host, band in (("lower",0.1,0.01),("upper",0.01,0.1)):
        out["sides"][side] = {}
        for name, d in (("near",.003),("far",.005)):
            layers = band_layers(d,host,band)
            out["sides"][side][name] = modal_admittances(layers,W,nmax).tolist()
    print(json.dumps(out, ensure_ascii=False, indent=2))

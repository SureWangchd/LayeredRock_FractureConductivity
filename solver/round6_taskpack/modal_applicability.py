"""Energy-based applicability diagnostic for an existing modal closure.

This does NOT solve COMSOL, replace the R5 solver, or certify continuum/mesh
errors.  It bounds a *specified finite modal correction* using a solved
reference field and known positive reference/target exterior admittances.

For two sides and n=0..N, psi_0=1 and psi_n=cos(n*pi*x/W).  N_0=W,
N_n=W/2. All lengths are metres and admittances S/m^2.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
import math
import numpy as np


@dataclass(frozen=True)
class ModalInterval:
    first_variation_S: float
    remainder_bound_S: float
    correction_lower_S: float
    correction_upper_S: float
    target_gain_lower_S: float
    target_gain_upper_S: float
    absolute_reference_error_bound_S: float
    relative_reference_error_bound: float | None
    absolute_mode_term_sum_S: float
    target_order: int

    def to_dict(self) -> dict:
        return asdict(self)


def modal_interval(
    reference_gain_S: float,
    reference_trace_V: np.ndarray,
    reference_Y_S_m2: np.ndarray,
    target_Y_S_m2: np.ndarray,
    *,
    width_m: float = 0.12,
    thickness_m: float = 1.0,
    voltage_V: float = 1.0,
    zero_mode_rtol: float = 1e-10,
) -> ModalInterval:
    """Bounds target minus reference gain within a common correction model.

    Arrays have shape (2,N+1); index zero is the constant mode and must match.
    `reference_trace_V` must come from the REFERENCE solve, not the target.
    The two models must share core, fixed electrode drive, true background
    conductance, FE space and high-mode tail. Differences outside the supplied
    modes are NOT covered. The caller must separately handle FE/mesh errors.

    Let Q=sum(t*Nn*dYn*v_ref,n^2)/DV^2 and
    B=sum(t*Nn*dYn^2/Y_target,n*v_ref,n^2)/DV^2. Then
         Q-B <= g_target-g_reference <= Q.
    Positivity of the full core/exterior energy, exact modal integrations and
    a sufficiently accurate stationary reference solve are required.
    """
    for name, value in [('width_m', width_m), ('thickness_m', thickness_m),
                        ('voltage_V', voltage_V)]:
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    if not math.isfinite(reference_gain_S):
        raise ValueError('reference_gain_S must be finite')
    v = np.asarray(reference_trace_V, dtype=float)
    yr = np.asarray(reference_Y_S_m2, dtype=float)
    yt = np.asarray(target_Y_S_m2, dtype=float)
    if v.ndim != 2 or v.shape[0] != 2 or v.shape[1] < 2:
        raise ValueError('Expected two sides with modes 0..N, N>=1')
    if v.shape != yr.shape or v.shape != yt.shape:
        raise ValueError('Trace/reference/target arrays must have equal shapes')
    if not (np.isfinite(v).all() and np.isfinite(yr).all() and np.isfinite(yt).all()):
        raise ValueError('Nonfinite array entries')
    if (yr <= 0).any() or (yt <= 0).any():
        raise ValueError('Complete reference and target admittances must be positive')
    if not np.allclose(yr[:, 0], yt[:, 0], rtol=zero_mode_rtol, atol=1e-13):
        raise ValueError('Changed zero modes/affine drive require a different formula')
    dy = yt[:, 1:] - yr[:, 1:]
    weight = thickness_m * (width_m / 2.0) * v[:, 1:]**2 / voltage_V**2
    terms = weight * dy
    Q = float(terms.sum())
    B = float((weight * dy**2 / yt[:, 1:]).sum())
    low, high = Q - B, Q
    lower, upper = reference_gain_S + low, reference_gain_S + high
    abs_bound = max(abs(low), abs(high))
    rel_bound = abs_bound / lower if lower > 0 else None
    return ModalInterval(Q, B, low, high, lower, upper, abs_bound,
                         rel_bound, float(np.abs(terms).sum()), v.shape[1] - 1)


def safe_for_relative_gain(interval: ModalInterval, tolerance: float = 0.02) -> bool:
    """Conservative model-level decision; 'False' means unconfirmed, not failed."""
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('tolerance must be positive')
    return interval.relative_reference_error_bound is not None and \
        interval.relative_reference_error_bound <= tolerance


def contrast_interval(
    lower: np.ndarray, upper: np.ndarray, coefficients: np.ndarray
) -> tuple[float, float]:
    """Signed linear contrast bounds. Shared errors are NOT assumed to cancel."""
    lo, hi, a = [np.asarray(x, dtype=float) for x in (lower, upper, coefficients)]
    if lo.shape != hi.shape or lo.shape != a.shape or not np.isfinite([lo, hi, a]).all():
        raise ValueError('Shapes must agree and inputs must be finite')
    if (lo > hi).any():
        raise ValueError('Lower limit exceeds upper limit')
    return (float(np.sum(np.where(a >= 0, a*lo, a*hi))),
            float(np.sum(np.where(a >= 0, a*hi, a*lo))))

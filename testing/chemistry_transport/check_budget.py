#!/usr/bin/env python3
"""Frozen double-precision acceptance checks for MICM gas transport.

Thresholds from docs/plans/pbl-transport.md ("Verification and completion
criteria"); they were frozen before any transport scheme was implemented, so
justify any later change in the plan. Shapes: m_d is (..., nz) dry layer mass
per area [kg m-2]; q is (..., nspecies, nz) [kg species / kg dry air].
"""
import numpy as np

MASS_RESID_PER_CALL = 1e-12      # isolated closed column, per call
MASS_RESID_1000_CALLS = 1e-10    # isolated closed column, 1,000 calls
FULL_RUN_RESID = 1e-8            # accumulated transport budget, full integration
UNIFORM_ATOL, UNIFORM_RTOL = 1e-20, 1e-12     # kg/kg; uniform-profile and permutation
PARALLEL_ATOL, PARALLEL_RTOL = 1e-18, 1e-10   # kg/kg; serial/8-rank, restart, short compiler


def column_mass(m_d, q):
    """Species mass per area [kg m-2], shape (..., nspecies)."""
    return np.sum(m_d[..., None, :] * q, axis=-1)


def relative_residual(m_d, q_before, q_after, boundary_flux=0.0, dt=0.0):
    """Mass change not explained by boundary_flux*dt [kg m-2 s-1 * s].

    Returns (relative, absolute); relative is normalized by
    max(initial species mass, sum(m_d)*1e-30), so a zero species stays finite.
    """
    before = column_mass(m_d, q_before)
    after = column_mass(m_d, q_after)
    norm = np.maximum(before, np.sum(m_d, axis=-1)[..., None] * 1e-30)
    absolute = np.abs(after - before - dt * boundary_flux)
    return absolute / norm, absolute


def nonnegative(q_after, q_initial):
    """True unless a value falls below -(1e-20 + 1e-12*max_column(q_initial))."""
    floor = -(1e-20 + 1e-12 * np.max(q_initial, axis=-1, keepdims=True))
    return bool(np.all(q_after >= floor))


def zero_species_stays_zero(q_after, q_initial):
    """Every species column that starts exactly zero stays exactly zero."""
    zero_column = np.all(q_initial == 0.0, axis=-1)
    return bool(np.all(q_after[zero_column] == 0.0))


def uniform_profile_unchanged(q_after, q_uniform):
    return bool(np.allclose(q_after, q_uniform, atol=UNIFORM_ATOL, rtol=UNIFORM_RTOL))


def disabled_matches_baseline(fields_new, fields_base):
    """Disabled transport reproduces every baseline field bit for bit."""
    for name, base in fields_base.items():
        if name not in fields_new:
            return False
        new = np.asarray(fields_new[name])
        base = np.asarray(base)
        if new.dtype != base.dtype or new.shape != base.shape or new.tobytes() != base.tobytes():
            return False
    return True

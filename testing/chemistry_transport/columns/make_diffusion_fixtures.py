#!/usr/bin/env python3
"""Generate reference fixtures for the closed-column implicit diffusion kernel.

The reference is independent of the Fortran Thomas sweep: it assembles the
flux-form backward-Euler system as a dense matrix and solves it exactly in
rational arithmetic (Gaussian elimination on fractions.Fraction), so the
stored solution is the exact solution rounded once to double,

    m_d(k) (q1(k) - q0(k)) / dt = G(k+1) (q1(k+1) - q1(k)) - G(k) (q1(k) - q1(k-1)),

where G(k) = rho*Kz/dz [kg m-2 s-1] at the interface below layer k, and
G(1) = G(nz+1) = 0 (closed column; MPAS level 1 at the ground).

Usage: make_diffusion_fixtures.py [OUTDIR]   (default: this directory)
Each case is written to diffusion_<name>.txt:
    # comment
    nz nt
    dt
    m_d(1:nz)
    G(1:nz+1)
    q0 rows (nt lines of nz values)
    q1 reference rows (nt lines)
"""
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


def solve(dt, m_d, G, q0):
    """Exact dense solve of the flux-form system for every species row of q0."""
    nz = m_d.size
    g = [Fraction(0)] + [Fraction(float(v)) for v in G[1:nz]] + [Fraction(0)]
    dt_f = Fraction(float(dt))
    m = [Fraction(float(v)) for v in m_d]
    out = np.empty_like(q0)
    for row in range(q0.shape[0]):
        A = [[Fraction(0)] * nz for _ in range(nz)]
        for k in range(nz):
            A[k][k] = m[k] / dt_f + g[k] + g[k + 1]
            if k > 0:
                A[k][k - 1] = -g[k]
            if k < nz - 1:
                A[k][k + 1] = -g[k + 1]
        b = [m[k] / dt_f * Fraction(float(q0[row, k])) for k in range(nz)]
        for col in range(nz):                                  # dense elimination, exact
            for r in range(col + 1, nz):
                if A[r][col] != 0:
                    f = A[r][col] / A[col][col]
                    for j in range(col, nz):
                        A[r][j] -= f * A[col][j]
                    b[r] -= f * b[col]
        x = [Fraction(0)] * nz
        for r in range(nz - 1, -1, -1):
            x[r] = (b[r] - sum(A[r][j] * x[j] for j in range(r + 1, nz))) / A[r][r]
        out[row] = [float(v) for v in x]
    return out


def cases():
    rng = np.random.default_rng(20261107)
    nz = 12
    m_d = 80.0 + 160.0 * rng.random(nz)                  # nonuniform dry layer mass [kg m-2]
    G = 0.05 + 2.0 * rng.random(nz + 1)                  # nonuniform exchange [kg m-2 s-1]
    G[0] = G[-1] = 0.0
    k = np.arange(nz)
    smooth = 4e-9 * (1.0 + 0.5 * np.sin(0.9 * k))
    step = np.where(k < nz // 2, 6e-9, 1e-12)
    yield "nonuniform", 30.0, m_d, G, np.vstack([smooth, np.zeros(nz), np.full(nz, 1e-30), step])

    pulses = np.zeros((3, nz))
    pulses[0, 0] = 5e-9
    pulses[1, nz // 2] = 5e-9
    pulses[2, nz - 1] = 5e-9
    yield "pulses", 60.0, m_d, G, pulses

    sharp = np.where(k < 3, 0.0, 2e-8)[None, :]
    yield "sharp_gradient", 10.0, m_d, G, sharp

    yield "zero_coefficients", 30.0, m_d, np.zeros(nz + 1), np.vstack([smooth, step])

    big = 40.0 + np.zeros(nz + 1)                        # dt*G/m_d = 960 with dt = 3600 s, m_d = 150
    big[0] = big[-1] = 0.0
    yield "large_timestep", 3600.0, m_d * 0 + 150.0, big, np.vstack([smooth, step, pulses[1]])

    yield "uniform", 30.0, m_d, G, np.full((2, nz), 3e-9)


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    for name, dt, m_d, G, q0 in cases():
        q1 = solve(dt, m_d, G, q0)
        lines = [f"# diffusion fixture '{name}', exact rational dense flux-form backward Euler",
                 f"{m_d.size} {q0.shape[0]}", f"{dt:.17e}",
                 " ".join(f"{v:.17e}" for v in m_d), " ".join(f"{v:.17e}" for v in G)]
        lines += [" ".join(f"{v:.17e}" for v in row) for row in q0]
        lines += [" ".join(f"{v:.17e}" for v in row) for row in q1]
        (outdir / f"diffusion_{name}.txt").write_text("\n".join(lines) + "\n")
        print(f"wrote diffusion_{name}.txt")


if __name__ == "__main__":
    main()

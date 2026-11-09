#!/usr/bin/env python3
"""In-model PBL gas mixing on the finer vertical grid (check_pbl_fine.sh; P6 of
docs/plans/pbl-transport.md).

  mixed ON_RUN OFF_RUN [--min-layers 8] [--band 2] [--dt 3] [--plumes] [--json PATH]
      At the final record:
      (1) the PBL spans at least MIN_LAYERS layers in the domain median: the
          number of layers whose centre lies below hpbl;
      (2) NO is mixed inside the PBL at the scheme's own rate: in each column,
          the normalized range (max - min)/mean of qNO over the layers below
          0.8*hpbl, as a domain median r (the initial NO falls off as
          exp(-z/300 m)). The mixing ln(r_off/r_on) must lie within a factor of
          BAND of ln(r_off/r_replay), where the replay diffuses the initial qNO
          offline with the scheme's own kzh from the output (each record's kzh
          held over the steps before it; closed column, implicit, dry-mass
          weighted, kzh on the bottom face of each layer as both YSU and MYNN
          write it). The replay must itself mix measurably (ln >= 0.1). It
          leaves out YSU's entrainment and MYNN's plumes, and for YSU it uses
          the heat coefficient, which the scalar one does not exceed inside the
          PBL, so it is a reference rate, not an exact one;
      (3) with --plumes, MYNN's updrafts are active: edmf_a > 0 somewhere in at
          least 1% of the columns at some record.
  --self-test
      Synthetic checks of the metrics.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def layers_below(zc, hpbl):
    """Number of layer centres below hpbl, per column: zc (cells, levels) [m], hpbl (cells) [m]."""
    return np.count_nonzero(zc < hpbl[:, None], axis=1)


def normalized_range(q, zc, hpbl, fraction=0.8):
    """(max - min)/mean of q over the layers below fraction*hpbl, per column (NaN with fewer than 2 layers)."""
    inside = zc < fraction * hpbl[:, None]
    out = np.full(q.shape[0], np.nan)
    for i in np.nonzero(inside.sum(axis=1) >= 2)[0]:
        v = q[i, inside[i]]
        out[i] = (v.max() - v.min()) / max(v.mean(), np.finfo(float).tiny)
    return out


def diffuse_step(q, kh, rho, dz, dt):
    """One implicit step of closed-column diffusion, per column: q, rho [kg m-3], dz [m] (cells, nz);
    kh [m2 s-1] (cells, nz) on the bottom face of each layer (kh[:, 0], the surface, is ignored).
    Interface exchange G = rhoz*kh/dzk between layers k-1 and k, rhoz the thickness-weighted density."""
    m = rho * dz
    dzk = 0.5 * (dz[:, 1:] + dz[:, :-1])
    rhoz = (rho[:, 1:] * dz[:, :-1] + rho[:, :-1] * dz[:, 1:]) / (dz[:, :-1] + dz[:, 1:])
    g = dt * rhoz * kh[:, 1:] / dzk                       # (cells, nz-1), interface below layer k = 1..nz-1
    flux = np.zeros((q.shape[0], q.shape[1] + 1))
    flux[:, 1:-1] = g * (q[:, :-1] - q[:, 1:])            # upward exchange through each interior interface
    rhs = flux[:, :-1] - flux[:, 1:]                      # dt * L(q): increment times m
    lower = np.zeros_like(q)
    upper = np.zeros_like(q)
    lower[:, 1:] = -g
    upper[:, :-1] = -g
    diag = m - lower - upper
    # (M - dt L) dq = dt L(q), Thomas algorithm along the column; K = 0 gives dq = 0 exactly
    cp = np.zeros_like(q)
    dp = np.zeros_like(q)
    cp[:, 0] = upper[:, 0] / diag[:, 0]
    dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for k in range(1, q.shape[1]):
        den = diag[:, k] - lower[:, k] * cp[:, k - 1]
        cp[:, k] = upper[:, k] / den
        dp[:, k] = (rhs[:, k] - lower[:, k] * dp[:, k - 1]) / den
    dq = np.zeros_like(q)
    dq[:, -1] = dp[:, -1]
    for k in range(q.shape[1] - 2, -1, -1):
        dq[:, k] = dp[:, k] - cp[:, k] * dq[:, k + 1]
    return q + dq


def record_steps(times, dt):
    """Model steps between consecutive records, from MPAS times 'YYYY-MM-DD_hh:mm:ss' within one month."""
    def seconds(s):
        date, clock = s.strip().split("_")
        h, mi, se = (int(x) for x in clock.split(":"))
        return int(date.split("-")[2]) * 86400 + h * 3600 + mi * 60 + se
    assert len({s.strip()[:7] for s in times}) == 1, "records span more than one month"
    t = [seconds(s) for s in times]
    steps = [(b - a) / dt for a, b in zip(t[:-1], t[1:])]
    assert all(s > 0 and abs(s - round(s)) < 1e-9 for s in steps), f"records not on the {dt}-s step: {times}"
    return [int(round(s)) for s in steps]


def mixing_agreement(r_off, r_on, r_replay, band, min_log=0.1):
    """(ok, fraction): fraction = ln(r_off/r_on) / ln(r_off/r_replay); ok when the replay mixes
    measurably (ln >= min_log) and the fraction lies in [1/band, band]."""
    l_rep = np.log(r_off / r_replay)
    frac = float(np.log(r_off / r_on) / l_rep) if l_rep > 0 else float("nan")
    return bool(l_rep >= min_log and 1.0 / band <= frac <= band), frac


def plume_fraction(edmf_a):
    """Largest fraction, over records, of columns with edmf_a > 0 at some level: edmf_a (time, cells, levels)."""
    return float(max((edmf_a[t] > 0.0).any(axis=1).mean() for t in range(edmf_a.shape[0])))


def cmd_mixed(args):
    from netCDF4 import Dataset
    failures, report = [], {}

    def read(run):
        with Dataset(Path(run) / "output.nc") as d:
            d.set_auto_mask(False)
            zg = np.asarray(d["zgrid"][:], dtype=float)
            zc = 0.5 * (zg[:, 1:] + zg[:, :-1]) - zg[:, :1]
            out = {"zc": zc, "hpbl": np.asarray(d["hpbl"][-1], dtype=float), "qNO": np.asarray(d["qNO"][-1], dtype=float)}
            if "edmf_a" in d.variables:
                out["edmf_a"] = np.asarray(d["edmf_a"][:], dtype=float)
            return out

    def replay(on_run, off_run, dt):
        """The off run's initial qNO diffused with the on run's kzh and rho, to the final record."""
        with Dataset(Path(on_run) / "output.nc") as d, Dataset(Path(off_run) / "output.nc") as o:
            d.set_auto_mask(False)
            o.set_auto_mask(False)
            dz = np.diff(np.asarray(d["zgrid"][:], dtype=float), axis=1)
            times = ["".join(c.decode() for c in row) for row in np.asarray(d["xtime"][:])]
            q = np.asarray(o["qNO"][0], dtype=float)
            for t, n in enumerate(record_steps(times, dt), start=1):
                kh, rho = np.asarray(d["kzh"][t], dtype=float), np.asarray(d["rho"][t], dtype=float)
                for _ in range(n):
                    q = diffuse_step(q, kh, rho, dz, dt)
            return q

    on, off = read(args.on), read(args.off)
    nlay = layers_below(on["zc"], on["hpbl"])
    report["pbl_layers_median"] = float(np.median(nlay))
    report["hpbl_median_m"] = float(np.median(on["hpbl"]))
    ok = report["pbl_layers_median"] >= args.min_layers
    print(f"{'PASS' if ok else 'FAIL'} the PBL spans {report['pbl_layers_median']:.0f} layers in the domain median "
          f"(hpbl {report['hpbl_median_m']:.0f} m; at least {args.min_layers} required)")
    if not ok:
        failures.append("pbl layers")
    r_on = np.nanmedian(normalized_range(on["qNO"], on["zc"], on["hpbl"]))
    r_off = np.nanmedian(normalized_range(off["qNO"], on["zc"], on["hpbl"]))
    q_rep = replay(args.on, args.off, args.dt)
    r_rep = np.nanmedian(normalized_range(q_rep, on["zc"], on["hpbl"]))
    ok, frac = mixing_agreement(r_off, r_on, r_rep, args.band)
    report.update(no_range_median_on=float(r_on), no_range_median_off=float(r_off),
                  no_range_median_replay=float(r_rep), mixing_fraction_of_replay=frac)
    print(f"{'PASS' if ok else 'FAIL'} NO inside the PBL: normalized range {r_on:.3g} with transport, {r_off:.3g} "
          f"without, {r_rep:.3g} in the local-K replay; mixing ln(off/on) {np.log(r_off / r_on):.3g} is "
          f"{frac:.3g} of the replay's {np.log(r_off / r_rep):.3g} ({1 / args.band:g} to {args.band:g} required)")
    if not ok:
        failures.append("mixing rate")
    if args.plumes:
        if "edmf_a" not in on:
            print("FAIL edmf_a is not in the output")
            failures.append("edmf_a missing")
        else:
            frac = plume_fraction(on["edmf_a"])
            report["plume_column_fraction"] = frac
            ok = frac >= 0.01
            print(f"{'PASS' if ok else 'FAIL'} MYNN updrafts active in {100 * frac:.2f}% of the columns "
                  f"(at least 1% required)")
            if not ok:
                failures.append("plumes")
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=1) + "\n")
    return 1 if failures else 0


def self_test():
    zc = np.tile(np.arange(25.0, 1500.0, 50.0), (3, 1))
    hpbl = np.array([500.0, 1000.0, 20.0])
    assert layers_below(zc, hpbl).tolist() == [10, 20, 0]
    q = np.exp(-zc / 300.0)
    mixed = np.full_like(q, 0.5)
    r = normalized_range(q, zc, hpbl)
    assert np.isnan(r[2]) and r[0] > 1.0, "a surface-confined profile has a large range; too shallow is NaN"
    assert np.allclose(normalized_range(mixed, zc, hpbl)[:2], 0.0), "a uniform profile has no range"
    a = np.zeros((2, 4, 3))
    a[1, 0, 1] = 0.02
    assert plume_fraction(a) == 0.25

    # replay kernel: closed column, implicit, K on the bottom face of each layer
    rng = np.random.default_rng(7)
    dz = rng.uniform(20.0, 80.0, (4, 30))
    rho = rng.uniform(0.9, 1.2, (4, 30))
    q0 = rng.uniform(0.0, 1.0, (4, 30))
    assert np.array_equal(diffuse_step(q0, np.zeros_like(q0), rho, dz, 3.0), q0), "K = 0 leaves q unchanged"
    q1 = diffuse_step(q0, rng.uniform(0.0, 200.0, (4, 30)), rho, dz, 3.0)
    m = rho * dz
    assert np.allclose((m * q1).sum(axis=1), (m * q0).sum(axis=1), rtol=1e-13, atol=0.0), "closed column conserves mass"
    # cos(pi z/H) on a uniform cell-centred grid with no-flux ends is an eigenvector of the discrete operator:
    # eigenvalue (4K/dz**2) sin**2(pi dz/(2H)); implicit Euler multiplies it by 1/(1 + dt*lambda) per step
    n, h, k, dt, steps = 100, 10.0, 50.0, 3.0, 200
    zc_u = (np.arange(n) + 0.5) * h
    mode = np.cos(np.pi * zc_u / (n * h))[None, :]
    q = mode.copy()
    for _ in range(steps):
        q = diffuse_step(q, np.full_like(q, k), np.ones_like(q), np.full_like(q, h), dt)
    face = np.zeros((1, 6))
    face[0, 3] = 100.0                                   # only the face between layers 2 and 3 exchanges
    qf = diffuse_step(np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]]), face, np.ones((1, 6)), np.full((1, 6), 10.0), 3.0)
    assert np.count_nonzero(qf[0, [0, 1, 4, 5]]) == 0 and qf[0, 2] > 0.0 and qf[0, 3] < 1.0, "kh[k] is the bottom face"
    lam = 4.0 * k / h**2 * np.sin(np.pi * h / (2.0 * n * h)) ** 2
    assert np.allclose(q, mode * (1.0 + dt * lam) ** (-steps), rtol=0.0, atol=1e-12), "discrete eigenmode decay"
    assert np.allclose(record_steps(["0000-01-01_18:00:00", "0000-01-01_18:01:00", "0000-01-02_00:00:30"], 3.0),
                       [20, 7190])

    # agreement of the in-model mixing with the replay: log-ratios within a factor of BAND
    assert mixing_agreement(2.0, 1.0, 1.0, 2.0)[0], "equal mixing agrees"
    assert not mixing_agreement(2.0, 2.0, 1.0, 2.0)[0], "no mixing in the model fails"
    assert not mixing_agreement(np.e**3, 1.0, np.e**2, 2.0)[0], "three times the replay's mixing fails"
    assert not mixing_agreement(1.05, 1.0, 1.0, 2.0)[0], "a replay without measurable mixing is no reference"
    ok, frac = mixing_agreement(np.e, np.e**0.4, 1.0, 2.0)
    assert ok and abs(frac - 0.6) < 1e-12
    print("check_pbl_fine self-test: PASS")


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if argv == ["--self-test"]:
        self_test()
        return 0
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("mixed")
    p.add_argument("on")
    p.add_argument("off")
    p.add_argument("--min-layers", type=int, default=8)
    p.add_argument("--band", type=float, default=2.0)
    p.add_argument("--dt", type=float, default=3.0)
    p.add_argument("--plumes", action="store_true")
    p.add_argument("--json")
    p.set_defaults(func=cmd_mixed)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

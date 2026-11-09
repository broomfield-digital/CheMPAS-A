#!/usr/bin/env python3
"""Self-test for check_transport_budget.py (plain asserts; run: python test_check_transport_budget.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_transport_budget as ctb  # noqa: E402


def fmt(x):
    """A value as the model writes it: es25.16e3, or a string such as gfortran's 'NaN' or 'Infinity'."""
    return x if isinstance(x, str) else f"{x:.16E}"


def line(t, gas, rel=1.0e-16, residual=1.0e-16, refresh=None, dt=3.0, columns=28080, qmin=0.0,
         before=1.0e6, after=1.0e6, increment=1.0e-10, qmax=1.0e-8):
    refresh = t if refresh is None else refresh
    return (f" [transport] budget PBL t={fmt(t)} gas={gas} scalar={10 + gas} before={fmt(before)} "
            f"after={fmt(after)} increment={fmt(increment)} rel={fmt(rel)} residual={fmt(residual)} "
            f"qmin={fmt(qmin)} qmax={fmt(qmax)} columns={columns} refresh={fmt(refresh)} dt={fmt(dt)}")


def log(steps=3, gases=3, interval=3.0, dt=3.0, **override):
    rows = []
    for s in range(steps):
        t = s * dt
        refresh = interval * (t // interval)
        for g in range(1, gases + 1):
            rows.append(line(t, g, refresh=refresh, dt=dt))
    return "\n".join(["some other log line"] + rows) + "\n"


def run(text, *args):
    return ctb.check(ctb.parse(text), ctb.build_parser().parse_args(["LOG", *args]))


def test_valid_log_passes():
    failures, summary = run(log(), "--expect-gases", "3", "--expect-steps", "3", "--dt", "3")
    assert failures == [], failures
    assert summary["steps"] == 3 and summary["lines"] == 9


def test_missing_gas_line_fails():
    text = "\n".join(l for l in log().splitlines() if "gas=2" not in l or "t=3.0" not in l)
    failures, _ = run(text, "--expect-gases", "3", "--expect-steps", "3")
    assert any("gas" in f for f in failures), failures


def test_relative_increment_over_threshold_fails():
    text = log() + line(9.0, 1, rel=2.0e-12) + "\n"
    failures, _ = run(text)
    assert any("rel" in f for f in failures), failures


def test_refresh_every_step_at_a_longer_interval_fails():
    failures, _ = run(log(steps=7, interval=3.0), "--interval", "9", "--dt", "3")
    assert any("refreshes" in f for f in failures), failures


def test_refresh_between_steps_fails():
    text = line(3.0, 1, refresh=1.5) + "\n"
    failures, _ = run(text, "--interval", "9", "--dt", "3")
    assert any("step time" in f for f in failures), failures


def test_interval_not_a_multiple_of_dt_passes():
    # 10-s interval, 3-s steps: the pbl alarm rings at the first step at or after 0, 10, 20, 30 s
    rows = []
    for s in range(14):
        t = 3.0 * s
        refresh = max(r for r in (0.0, 12.0, 21.0, 30.0) if r <= t)
        rows.append(line(t, 1, refresh=refresh))
    failures, _ = run("\n".join(rows) + "\n", "--interval", "10", "--dt", "3", "--expect-steps", "14")
    assert failures == [], failures


def test_operator_held_for_a_whole_interval_fails():
    text = line(9.0, 1, refresh=0.0) + "\n"
    failures, _ = run(text, "--interval", "9")
    assert any("older" in f for f in failures), failures


def test_held_operator_within_its_interval_passes():
    failures, _ = run(log(steps=7, interval=9.0), "--interval", "9", "--expect-steps", "7", "--expect-gases", "3")
    assert failures == [], failures


def test_wrong_dt_fails():
    failures, _ = run(log(dt=3.0), "--dt", "6")
    assert any("dt" in f for f in failures), failures


def test_nan_budget_fails():
    # gfortran's MAX/MIN drop NaN, so the reduced residual and qmin can look clean while the masses are NaN
    text = log() + line(6.0, 1, after="NaN", increment="NaN", rel="NaN", residual="NaN", qmin="NaN", refresh=6.0) + "\n"
    failures, _ = run(text, "--dt", "3", "--interval", "3", "--min-q", "0")
    assert any("not finite" in f for f in failures), failures


def test_nan_mass_with_clean_rel_fails():
    text = log() + line(6.0, 1, before="NaN", refresh=6.0) + "\n"
    failures, _ = run(text)
    assert any("before" in f and "not finite" in f for f in failures), failures


def test_infinite_qmax_fails():
    text = log() + line(6.0, 1, qmax="Infinity", refresh=6.0) + "\n"
    failures, _ = run(text)
    assert any("qmax" in f and "not finite" in f for f in failures), failures


def test_increment_inconsistent_with_masses_fails():
    # the increment is summed column by column, independently of the two totals
    text = log() + line(6.0, 1, before=1.0e6, after=1.0e6 + 1.0e-3, increment=0.0, rel=0.0, refresh=6.0) + "\n"
    failures, _ = run(text)
    assert any("after - before" in f for f in failures), failures


def test_absent_gas_passes():
    # before lightning starts, NO is zero everywhere: before = after = increment = 0
    text = log() + line(6.0, 1, before=0.0, after=0.0, increment=0.0, rel=0.0, qmax=0.0, refresh=6.0) + "\n"
    failures, _ = run(text)
    assert failures == [], failures


def test_no_lines_fails():
    failures, _ = run("nothing here\n")
    assert failures, "an empty log must fail"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"test_check_transport_budget: {len(tests)} passed")

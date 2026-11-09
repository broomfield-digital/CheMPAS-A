#!/usr/bin/env python3
"""Self-test for check_transport_budget.py (plain asserts; run: python test_check_transport_budget.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_transport_budget as ctb  # noqa: E402


def line(t, gas, rel=1.0e-16, residual=1.0e-16, refresh=None, dt=3.0, columns=28080, qmin=0.0):
    refresh = t if refresh is None else refresh
    return (f" [transport] budget PBL t={t:.16E} gas={gas} scalar={10 + gas} before=1.0E+06 after=1.0E+06 "
            f"increment=1.0E-10 rel={rel:.16E} residual={residual:.16E} qmin={qmin:.16E} qmax=1.0E-08 "
            f"columns={columns} refresh={refresh:.16E} dt={dt:.16E}")


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


def test_no_lines_fails():
    failures, _ = run("nothing here\n")
    assert failures, "an empty log must fail"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"test_check_transport_budget: {len(tests)} passed")

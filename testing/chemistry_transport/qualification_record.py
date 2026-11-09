#!/usr/bin/env python3
"""Provenance, summary and manifest of a PBL gas transport qualification run
(run_pbl_qualification.sh; P6 of docs/plans/pbl-transport.md).

  start RUN_ROOT REPO MODEL BASELINE [--reference NAME=PATH ...] [--source-model PATH]
        [--fine-grid DIR]
      Before the first check: write RUN_ROOT/provenance-start.json with the
      source revision, branch and tracked changes, the SHA-256 of the model and
      of each reference executable (the copies the checks run), the executable
      they were copied from (recorded, not checked for drift: a rebuild there
      cannot reach the copies), the toolchain, the baseline and its hashes, and
      the fine grid's inputs.
  finish RUN_ROOT REPO
      After the last check: read RUN_ROOT/results.tsv (check, PASS or FAIL,
      seconds, log) and each check's log, record the source and the executables
      again, and write RUN_ROOT/summary.json and RUN_ROOT/manifest.json. The
      status is "fail" if a check failed or the source or an executable changed
      during the run, "incomplete" if a check skipped part of itself (a line
      starting with SKIP in its log), and "pass" otherwise. The exit status is 0
      only for "pass".
  --self-test
      Synthetic checks of the summary and the drift detection.
"""
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sh(*cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def first_line(text):
    return text.splitlines()[0] if text else None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_state(repo):
    """What the run must not change: revision, branch and tracked changes of the checkout."""
    return {
        "revision": sh("git", "-C", str(repo), "rev-parse", "HEAD"),
        "branch": sh("git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"),
        "tracked_changes": sh("git", "-C", str(repo), "status", "--short", "--untracked-files=no") or "",
    }


def executables_state(executables):
    """{name: {path, sha256}} for {name: path}."""
    return {name: {"path": str(path), "sha256": sha256(path)} for name, path in executables.items()}


def drift(start, end):
    """Differences between the source and executable states at the start and the end, as messages."""
    out = [f"source {k} changed during the run: {start['source'][k]!r} -> {end['source'][k]!r}"
           for k in start["source"] if start["source"][k] != end["source"].get(k)]
    for name, e in start["executables"].items():
        now = end["executables"].get(name)
        if now is None or now["sha256"] != e["sha256"]:
            out.append(f"executable {name} ({e['path']}) changed during the run")
    return out


def skipped_lines(log_text):
    return [line.strip() for line in log_text.splitlines() if line.startswith("SKIP")]


def summarize(rows, logs, drift_messages):
    """Summary of rows [(check, result, seconds, log path)] with logs {log path: text}."""
    checks = []
    for name, result, seconds, log in rows:
        checks.append({"check": name, "result": result, "elapsed_s": int(seconds), "log": log,
                       "skipped": skipped_lines(logs.get(log, ""))})
    if drift_messages or any(c["result"] != "PASS" for c in checks) or not checks:
        status = "fail"
    elif any(c["skipped"] for c in checks):
        status = "incomplete"
    else:
        status = "pass"
    return {"status": status, "provenance_drift": drift_messages, "checks": checks}


def cmd_start(args):
    run_root, repo, model, baseline = map(Path, args[:4])
    references, fine, source_model = {}, None, None
    rest = args[4:]
    while rest:
        flag, value, rest = rest[0], rest[1], rest[2:]
        if flag == "--reference":
            name, _, path = value.partition("=")
            references[name] = Path(path)
        elif flag == "--fine-grid":
            fine = Path(value)
        elif flag == "--source-model":
            source_model = {"path": value, "sha256": sha256(value)}
        else:
            raise SystemExit(f"unknown option {flag}")
    base_manifest = json.loads((baseline / "manifest.json").read_text()) if (baseline / "manifest.json").exists() else {}
    fine_inputs = {}
    if fine is not None:
        for name in ("convective_init.nc", "stable_init.nc", "supercell_zeta_levels.txt"):
            if (fine / name).exists():
                fine_inputs[name] = sha256(fine / name)
        rev = fine / "provenance" / "chempas_revision"
        fine_inputs["prepared_at_revision"] = rev.read_text().strip() if rev.exists() else None
    record = {
        "date_utc_start": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": platform.node(),
        "run_root": str(run_root),
        "source": source_state(repo),
        "executables": executables_state({"model": model, **references}),
        "model_copied_from": source_model,
        "toolchain": {
            "platform": f"{platform.system()} {platform.machine()}",
            "compiler": first_line(sh("gfortran", "--version")),
            "mpi": first_line(sh("mpirun", "--version")),
            "netcdf_c": sh("nc-config", "--version"),
            "netcdf_fortran": sh("nf-config", "--version"),
            "musica_prefix": os.environ.get("MUSICA_PREFIX"),
        },
        "baseline": {
            "directory": str(baseline),
            "supercell_init_sha256": sha256(baseline / "supercell_init.nc"),
            "chemistry": base_manifest.get("chemistry", {}),
            "toolchain": base_manifest.get("toolchain", {}),
        },
        "fine_grid": {"directory": str(fine) if fine else None, "inputs_sha256": fine_inputs},
    }
    (run_root / "provenance-start.json").write_text(json.dumps(record, indent=1) + "\n")
    return 0


def cmd_finish(args):
    run_root, repo = map(Path, args[:2])
    start = json.loads((run_root / "provenance-start.json").read_text())
    end = {"source": source_state(repo),
           "executables": {name: {"path": e["path"], "sha256": sha256(e["path"]) if Path(e["path"]).exists() else None}
                           for name, e in start["executables"].items()}}
    rows = [line.split("\t") for line in (run_root / "results.tsv").read_text().splitlines() if line]
    logs = {log: Path(log).read_text(errors="replace") for _, _, _, log in rows if Path(log).exists()}
    summary = summarize(rows, logs, drift(start, end))
    (run_root / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    manifest = {
        "status": summary["status"],
        "case": "MICM gas PBL transport qualification (docs/plans/pbl-transport.md, P5-P6)",
        **start,
        "date_utc_end": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "provenance_drift": summary["provenance_drift"],
        "checks": {c["check"]: c["result"] + (" (skipped part)" if c["skipped"] else "") for c in summary["checks"]},
    }
    (run_root / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    for message in summary["provenance_drift"]:
        print(f"FAIL {message}")
    for c in summary["checks"]:
        for line in c["skipped"]:
            print(f"{c['check']}: {line}")
    passed = sum(c["result"] == "PASS" for c in summary["checks"])
    print(f"{summary['status'].upper()}: {passed} of {len(summary['checks'])} checks passed")
    return 0 if summary["status"] == "pass" else 1


def self_test():
    rows = [("a", "PASS", "3", "a.log"), ("b", "PASS", "5", "b.log")]
    clean = {"a.log": "PASS one\n", "b.log": "PASS two\n"}
    assert summarize(rows, clean, [])["status"] == "pass"
    skipped = {"a.log": "PASS one\n", "b.log": "PASS two\nSKIP no PRE_P5_MODEL given\n"}
    s = summarize(rows, skipped, [])
    assert s["status"] == "incomplete", "a check that skipped part of itself is not a pass"
    assert s["checks"][1]["skipped"] == ["SKIP no PRE_P5_MODEL given"]
    assert summarize([("a", "FAIL", "3", "a.log")] + rows[1:], skipped, [])["status"] == "fail"
    assert summarize(rows, clean, ["executable model changed during the run"])["status"] == "fail"
    assert summarize([], {}, [])["status"] == "fail", "no checks is not a pass"

    start = {"source": {"revision": "abc", "branch": "b", "tracked_changes": ""},
             "executables": {"model": {"path": "/m", "sha256": "1"}, "pre_p5": {"path": "/p", "sha256": "2"}}}
    same = json.loads(json.dumps(start))
    assert drift(start, same) == []
    moved = json.loads(json.dumps(start))
    moved["source"]["revision"] = "def"
    moved["source"]["tracked_changes"] = " M src/x.F"
    moved["executables"]["model"]["sha256"] = "9"
    d = drift(start, moved)
    assert len(d) == 3 and any("revision" in m for m in d) and any("model" in m for m in d), d
    gone = json.loads(json.dumps(start))
    del gone["executables"]["pre_p5"]
    assert any("pre_p5" in m for m in drift(start, gone))
    print("qualification_record self-test: PASS")


def main(argv):
    if argv == ["--self-test"]:
        self_test()
        return 0
    if len(argv) >= 5 and argv[0] == "start":
        return cmd_start(argv[1:])
    if len(argv) == 3 and argv[0] == "finish":
        return cmd_finish(argv[1:])
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

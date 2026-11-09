#!/usr/bin/env bash
# Full qualification of MICM gas PBL transport with the GNU toolchain on eight
# ranks (P6 of docs/plans/pbl-transport.md): every unit test and in-model check,
# in order, with a summary and a manifest of the build and environment.
#
# Usage: run_pbl_qualification.sh BASELINE_DIR RUN_ROOT [MODEL] [--dry-run]
#
# Each check writes its runs under RUN_ROOT/<check> and its log to
# RUN_ROOT/<check>.log. --dry-run lists the commands and checks that every
# script exists. Afterwards RUN_ROOT holds:
#   summary.json    each check: PASS or FAIL, elapsed seconds, log
#   manifest.json   revision, source status, executable SHA-256, compiler,
#                   MPI, NetCDF and MUSICA versions, baseline and its hashes
# The 60-s smoke run is compared with the tracked P1 reference
# (results/supercell-lnox-gnu-20261107/smoke-60s-summary.json). The fine grid
# is prepared under RUN_ROOT/fine-grid unless CHEMPAS_FINE_GRID names an
# existing one (prepare_fine_grid.sh). RUN_ROOT must not exist.
set -euo pipefail
dry=0
args=()
for a in "$@"; do
    if [[ $a == --dry-run ]]; then dry=1; else args+=("$a"); fi
done
if (( ${#args[@]} < 2 || ${#args[@]} > 3 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL] [--dry-run]" >&2
    exit 2
fi
baseline=$(cd "${args[0]}" && pwd)
run_root=${args[1]}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${args[2]:-$repo/atmosphere_model}")
test -x "$model"
if (( ! dry )) && [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
fine=${CHEMPAS_FINE_GRID:-$run_root/fine-grid}
smoke_reference=$here/results/supercell-lnox-gnu-20261107/smoke-60s-summary.json

# name | command (run from the repository root; $B baseline, $R run root, $M model)
checks=(
    "python-unit|for t in testing/chemistry_transport/test_check_*.py; do python \"\$t\" || exit 1; done && python testing/chemistry_transport/check_pbl_reactive.py --self-test && python testing/chemistry_transport/check_pbl_fine.py --self-test && python testing/chemistry_transport/make_fine_zeta_levels.py --self-test"
    "fortran-unit|bash testing/chemistry_transport/unit/run_unit_tests.sh"
    "mmm-patches|bash testing/chemistry_transport/test_apply_mmm_patches.sh"
    "wiring|testing/chemistry_transport/check_transport_wiring.sh \"\$B\" \"\$R/wiring\" \"\$M\""
    "realized-increment|testing/chemistry_transport/check_realized_increment.sh \"\$B\" \"\$R/realized-increment\""
    "ysu|testing/chemistry_transport/check_pbl_ysu.sh \"\$B\" \"\$R/ysu\" \"\$M\""
    "mynn|testing/chemistry_transport/check_pbl_mynn.sh \"\$B\" \"\$R/mynn\" \"\$M\""
    "cadence|testing/chemistry_transport/check_pbl_cadence.sh \"\$B\" \"\$R/cadence\" \"\$M\""
    "intervals|testing/chemistry_transport/check_pbl_intervals.sh \"\$B\" \"\$R/intervals\" \"\$M\""
    "restart|testing/chemistry_transport/check_pbl_restart.sh \"\$B\" \"\$R/restart\" \"\$M\""
    "chemistry-intervals|testing/chemistry_transport/check_pbl_chemistry_intervals.sh \"\$B\" \"\$R/chemistry-intervals\" \"\$M\""
    "mpi|testing/chemistry_transport/check_pbl_mpi.sh \"\$B\" \"\$R/mpi\" \"\$M\""
    "reactive|testing/chemistry_transport/check_pbl_reactive.sh \"\$B\" \"\$R/reactive\" \"\$M\""
    "fine-grid|{ [[ -s \"\$F/convective_init.nc\" ]] || testing/chemistry_transport/prepare_fine_grid.sh \"\$B\" \"\$F\"; } && testing/chemistry_transport/check_pbl_fine.sh \"\$B\" \"\$F\" \"\$R/fine\" \"\$M\""
    "smoke-60s|CHEMPAS_SMOKE_REFERENCE=\"\$S\" bash testing/chemistry_transport/run_supercell_smoke.sh \"\$B\" \"\$R/smoke-60s\" 60 \"\$M\""
)

if (( dry )); then
    missing=0
    for entry in "${checks[@]}"; do
        name=${entry%%|*}
        cmd=${entry#*|}
        echo "$name: $cmd"
        for script in $(grep -oE 'testing/chemistry_transport/[A-Za-z0-9_./-]+\.(sh|py)' <<< "$cmd" | sort -u); do
            if [[ $script == *'*'* ]]; then continue; fi
            [[ -f "$repo/$script" ]] || { echo "  MISSING $script"; missing=1; }
        done
    done
    [[ -f "$smoke_reference" ]] || { echo "  MISSING $smoke_reference"; missing=1; }
    (( missing == 0 )) && echo "dry run: every script exists"
    exit "$missing"
fi

mkdir -p "$run_root"
run_root=$(cd "$run_root" && pwd)
[[ $fine == /* ]] || fine=$run_root/fine-grid
source "$baseline/environment.sh"
export B=$baseline R=$run_root M=$model F=$fine S=$smoke_reference
results=$run_root/results.tsv
: > "$results"
cd "$repo"
for entry in "${checks[@]}"; do
    name=${entry%%|*}
    cmd=${entry#*|}
    echo "== $name"
    start=$(date +%s)
    set +e
    bash -c "set -euo pipefail; $cmd" > "$run_root/$name.log" 2>&1
    rc=$?
    set -e
    elapsed=$(( $(date +%s) - start ))
    result=PASS
    (( rc == 0 )) || result=FAIL
    printf '%s\t%s\t%s\t%s\n' "$name" "$result" "$elapsed" "$run_root/$name.log" >> "$results"
    echo "$result $name (${elapsed} s); log: $run_root/$name.log"
done

python3 - "$run_root" "$repo" "$model" "$baseline" <<'PY'
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

run_root, repo, model, baseline = map(Path, sys.argv[1:])


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


rows = [line.split("\t") for line in (run_root / "results.tsv").read_text().splitlines() if line]
checks = [{"check": n, "result": r, "elapsed_s": int(e), "log": log} for n, r, e, log in rows]
summary = {"status": "pass" if all(c["result"] == "PASS" for c in checks) else "fail", "checks": checks}
(run_root / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")

base_manifest = json.loads((baseline / "manifest.json").read_text()) if (baseline / "manifest.json").exists() else {}
manifest = {
    "status": summary["status"],
    "case": "MICM gas PBL transport qualification (docs/plans/pbl-transport.md, P5-P6)",
    "date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "host": platform.node(),
    "run_root": str(run_root),
    "source": {
        "revision": sh("git", "-C", str(repo), "rev-parse", "HEAD"),
        "branch": sh("git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"),
        "tracked_changes": sh("git", "-C", str(repo), "status", "--short", "--untracked-files=no") or "",
    },
    "executable": {"path": str(model), "sha256": sha256(model)},
    "toolchain": {
        "platform": f"{platform.system()} {platform.machine()}",
        "compiler": first_line(sh("gfortran", "--version")),
        "mpi": first_line(sh("mpirun", "--version")),
        "netcdf_c": sh("nc-config", "--version"),
        "netcdf_fortran": sh("nf-config", "--version"),
        "musica_prefix": os.environ.get("MUSICA_PREFIX"),
        "precision": "double",
        "openmp": False,
        "ranks": 8,
    },
    "baseline": {
        "directory": str(baseline),
        "supercell_init_sha256": sha256(baseline / "supercell_init.nc"),
        "chemistry": base_manifest.get("chemistry", {}),
        "toolchain": base_manifest.get("toolchain", {}),
    },
    "checks": {c["check"]: c["result"] for c in checks},
}
(run_root / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
print(f"{summary['status'].upper()}: {sum(c['result'] == 'PASS' for c in checks)} of {len(checks)} checks passed")
PY
grep -q $'\tFAIL\t' "$results" && exit 1
exit 0

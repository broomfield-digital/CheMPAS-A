#!/usr/bin/env bash
# Reuse frozen inputs, but test the current checkout's freshly built executable.
# Usage: run_supercell_smoke.sh BASELINE_DIR NEW_RUN_DIR [60|600] [MODEL]
# Optional CHEMPAS_SMOKE_REFERENCE supplies a same-toolchain JSON reference.
set -euo pipefail
if (( $# < 2 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR NEW_RUN_DIR [60|600] [MODEL]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_dir=$2
duration=${3:-60}
script_dir=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$script_dir/../.." && pwd)
model=${4:-$repo/atmosphere_model}
reference=${CHEMPAS_SMOKE_REFERENCE:-}
if [[ -n "$reference" ]]; then
    reference=$(realpath "$reference")
    test -f "$reference"
fi
case "$duration" in 60|600) ;; *) echo "Smoke duration must be 60 or 600 seconds" >&2; exit 2 ;; esac
test -x "$model"
test -f "$baseline/supercell_init.nc"
test -f "$baseline/environment.sh"
if [[ -e "$run_dir" ]]; then
    echo "Use a new run directory; existing runs are preserved: $run_dir" >&2
    exit 2
fi
mkdir -p "$run_dir/bin" "$run_dir/provenance"
run_dir=$(cd "$run_dir" && pwd)
for name in environment.sh namelist.atmosphere streams.atmosphere stream_list.atmosphere.output \
            supercell.graph.info.part.8 lnox_o3.yaml tuvx_no2.json tuvx_upper_atm.csv; do
    cp "$baseline/$name" "$run_dir/$name"
done
cp --reflink=auto "$baseline/supercell_init.nc" "$run_dir/"
cp "$model" "$run_dir/bin/atmosphere_model"
ln -s bin/atmosphere_model "$run_dir/atmosphere_model"
ln -s "$baseline/data" "$run_dir/data"
shopt -s nullglob
for path in "$baseline"/*.TBL "$baseline"/*DATA*; do
    ln -s "$path" "$run_dir/$(basename "$path")"
done
source "$run_dir/environment.sh"
python - "$run_dir" "$duration" <<'PY'
from pathlib import Path
import re
import sys
run, seconds = Path(sys.argv[1]), int(sys.argv[2])
duration = f"00:{seconds // 60:02d}:{seconds % 60:02d}"
path = run / 'namelist.atmosphere'
text, count = re.subn(r"config_run_duration\s*=\s*'[^']*'", f"config_run_duration = '{duration}'", path.read_text())
if count != 1:
    raise SystemExit('Expected exactly one run duration')
path.write_text(text)
path = run / 'streams.atmosphere'
interval = seconds // 2
path.write_text(path.read_text().replace('output_interval="00:10:00"', f'output_interval="00:{interval // 60:02d}:{interval % 60:02d}"'))
PY
git -C "$repo" rev-parse HEAD > "$run_dir/provenance/chempas_revision"
git -C "$repo" status --short > "$run_dir/provenance/source-status"
git -C "$repo" diff --binary -- src Makefile > "$run_dir/provenance/source.patch"
cd "$run_dir"
sha256sum bin/atmosphere_model supercell_init.nc namelist.atmosphere streams.atmosphere \
    lnox_o3.yaml tuvx_no2.json tuvx_upper_atm.csv > provenance/smoke-input-sha256.txt
date -u +%FT%TZ > provenance/run-started-utc
set +e
/usr/bin/time -v -o provenance/run-time.txt mpiexec -n 8 ./atmosphere_model > run.out 2>&1
run_status=$?
set -e
printf '%s\n' "$run_status" > provenance/run-exit-status
date -u +%FT%TZ > provenance/run-finished-utc
if (( run_status != 0 )); then
    echo "MPI run failed ($run_status); see $run_dir/run.out" >&2
    exit "$run_status"
fi
check_args=("$run_dir" --duration-seconds "$duration")
if [[ -n "$reference" ]]; then
    check_args+=(--reference "$reference")
fi
python "$script_dir/check_supercell_baseline.py" "${check_args[@]}"

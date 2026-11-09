#!/usr/bin/env bash
# Full qualification of MICM gas PBL transport with the GNU toolchain on eight
# ranks (P6 of docs/plans/pbl-transport.md): every unit test and in-model check,
# in order, with a summary and a manifest of the build and environment.
#
# Usage: run_pbl_qualification.sh BASELINE_DIR RUN_ROOT [MODEL] [--dry-run]
#
# Each check writes its runs under RUN_ROOT/<check> and its log to
# RUN_ROOT/<check>.log. --dry-run lists the commands and checks that every
# script and reference exists.
#
# The model and the reference executables are copied into RUN_ROOT before the
# first check, and every check runs the copies. The base comparisons need the
# executables before each stage, named by environment variables:
#   CHEMPAS_PRE_P3_MODEL  before P3 (check_pbl_ysu.sh: transport off, bitwise)
#   CHEMPAS_PRE_P4_MODEL  before P4, with the MYNN snow-content fix
#                         (check_pbl_mynn.sh)
#   CHEMPAS_PRE_P5_MODEL  before P5 (check_pbl_cadence.sh: no interval, bitwise)
# Without one, that check skips its comparison and the run is "incomplete".
#
# qualification_record.py records the source, executables, toolchain, baseline
# and fine grid at the start, and again at the end. Afterwards RUN_ROOT holds:
#   summary.json    the status (pass, incomplete or fail) and each check:
#                   PASS or FAIL, elapsed seconds, log, any SKIP lines
#   manifest.json   the start record, the end time, any change to the source or
#                   an executable during the run (a failure), and the results
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

declare -A references=(
    [pre_p3]=${CHEMPAS_PRE_P3_MODEL:-}
    [pre_p4]=${CHEMPAS_PRE_P4_MODEL:-}
    [pre_p5]=${CHEMPAS_PRE_P5_MODEL:-}
)

# name | command (run from the repository root; $B baseline, $R run root, $M model,
# $P3, $P4, $P5 the reference executables or empty)
checks=(
    "python-unit|for t in testing/chemistry_transport/test_check_*.py; do python \"\$t\" || exit 1; done && python testing/chemistry_transport/check_pbl_reactive.py --self-test && python testing/chemistry_transport/check_pbl_fine.py --self-test && python testing/chemistry_transport/make_fine_zeta_levels.py --self-test && python testing/chemistry_transport/qualification_record.py --self-test"
    "fortran-unit|bash testing/chemistry_transport/unit/run_unit_tests.sh"
    "mmm-patches|bash testing/chemistry_transport/test_apply_mmm_patches.sh"
    "wiring|testing/chemistry_transport/check_transport_wiring.sh \"\$B\" \"\$R/wiring\" \"\$M\""
    "realized-increment|testing/chemistry_transport/check_realized_increment.sh \"\$B\" \"\$R/realized-increment\""
    "ysu|testing/chemistry_transport/check_pbl_ysu.sh \"\$B\" \"\$R/ysu\" \"\$M\" \${P3:+\"\$P3\"}"
    "mynn|testing/chemistry_transport/check_pbl_mynn.sh \"\$B\" \"\$R/mynn\" \"\$M\" \${P4:+\"\$P4\"}"
    "cadence|testing/chemistry_transport/check_pbl_cadence.sh \"\$B\" \"\$R/cadence\" \"\$M\" \${P5:+\"\$P5\"}"
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
    for name in pre_p3 pre_p4 pre_p5; do
        path=${references[$name]}
        if [[ -z "$path" ]]; then
            echo "  NO REFERENCE $name: its base comparison will be skipped (the run is then incomplete)"
        elif [[ ! -x "$path" ]]; then
            echo "  MISSING reference $name: $path"; missing=1
        else
            echo "  reference $name: $path"
        fi
    done
    (( missing == 0 )) && echo "dry run: every script exists"
    exit "$missing"
fi

for name in pre_p3 pre_p4 pre_p5; do
    [[ -z "${references[$name]}" || -x "${references[$name]}" ]] || { echo "Not executable: ${references[$name]}" >&2; exit 2; }
done
mkdir -p "$run_root/executables"
run_root=$(cd "$run_root" && pwd)
[[ $fine == /* ]] || fine=$run_root/fine-grid
source "$baseline/environment.sh"
# the checks run copies, so a rebuild in the checkout cannot reach them mid-run
cp -p "$model" "$run_root/executables/atmosphere_model"
record=(--source-model "$model")
P3= P4= P5=
for name in pre_p3 pre_p4 pre_p5; do
    [[ -n "${references[$name]}" ]] || continue
    cp -p "${references[$name]}" "$run_root/executables/$name-atmosphere_model"
    record+=(--reference "$name=$run_root/executables/$name-atmosphere_model")
done
[[ -n "${references[pre_p3]}" ]] && P3=$run_root/executables/pre_p3-atmosphere_model
[[ -n "${references[pre_p4]}" ]] && P4=$run_root/executables/pre_p4-atmosphere_model
[[ -n "${references[pre_p5]}" ]] && P5=$run_root/executables/pre_p5-atmosphere_model
model=$run_root/executables/atmosphere_model
[[ -d "$fine" ]] && record+=(--fine-grid "$fine")
python3 "$here/qualification_record.py" start "$run_root" "$repo" "$model" "$baseline" "${record[@]}"
export B=$baseline R=$run_root M=$model F=$fine S=$smoke_reference P3 P4 P5
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

set +e
python3 "$here/qualification_record.py" finish "$run_root" "$repo"
rc=$?
set -e
exit "$rc"

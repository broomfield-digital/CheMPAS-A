#!/usr/bin/env bash
# In-model qualification of YSU PBL mixing of MICM gases (P3 of
# docs/plans/pbl-transport.md), on the supercell with YSU and the
# Monin-Obukhov surface layer, an inert mechanism (titration A = 0, j = 0,
# TUV-x off) and lightning off. Gas profiles with vertical structure come from
# make_pbl_init.py; the gases are passive.
#
# Usage: check_pbl_ysu.sh BASELINE_DIR RUN_ROOT [MODEL [PRE_P3_MODEL]]
#
# Cases:
#   gate-mynn, gate-interval   transport requested with MYNN or a nondefault
#                              PBL interval stops the model with the gate's
#                              message
#   gate-no-sfclayer           YSU without a surface layer stops in MPAS's own
#                              physics_namelist_check, before transport init
#                              (the gas gate's message is unit-tested)
#   gate-pbl-off               transport requested with PBL off runs as a
#                              logged zero-transport case
#   stable-on / stable-off     baseline surface (288 K, colder than the air):
#                              with and without gas transport, 300 s
#   convective-on / -off       skintemp = sst = 305 K: upward heat flux
#   convective-topdown-on/-off the convective case with config_ysu_pblmix =
#                              true. That mode raises YSU's PBL height with the
#                              liquid-ice potential temperature, changing the
#                              meteorology; its radiatively driven cloud-top
#                              entrainment needs cooling, and the supercell has
#                              no radiation (the column reference tests cover
#                              that branch)
#   stable-off-pre-p3          stable-off with PRE_P3_MODEL (the executable
#                              before P3): all output bitwise identical
# Checks are in check_pbl_ysu.py: each on/off pair has bitwise identical
# non-gas fields, hpbl > 0, conserved gas mass and down-gradient mixing of
# every gas; the convective case mixes the gases whose gradients lie inside
# the boundary layer (NO, O3) more strongly than the stable case. RUN_ROOT must
# not exist.
set -euo pipefail
if (( $# < 2 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR RUN_ROOT [MODEL [PRE_P3_MODEL]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_root=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${3:-$repo/atmosphere_model}")
pre_p3_model=${4:+$(realpath "$4")}
test -x "$model"
[[ -z "$pre_p3_model" ]] || test -x "$pre_p3_model"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" ${pre_p3_model:+"$pre_p3_model"} > "$run_root/provenance/model-sha256.txt"

python "$here/make_pbl_init.py" "$baseline/supercell_init.nc" "$run_root/stable_init.nc"
python "$here/make_pbl_init.py" "$baseline/supercell_init.nc" "$run_root/convective_init.nc" --skin-temperature 305
python - "$baseline/lnox_o3.yaml" "$run_root/inert.yaml" <<'PY'
import re
import sys
from pathlib import Path
text = Path(sys.argv[1]).read_text()
text, count = re.subn(r'^(\s*A:\s*)\S+', r'\g<1>0.0', text, flags=re.M)
assert count == 1, 'expected one Arrhenius A'
Path(sys.argv[2]).write_text(text)
PY

stage() {
    # stage CASE INIT MODEL DURATION KEY=VALUE...  (namelist options; transport=default leaves it unset)
    local case_dir=$run_root/$1 init=$2 exe=$3 duration=$4
    mkdir -p "$case_dir"
    for name in namelist.atmosphere streams.atmosphere stream_list.atmosphere.output \
                supercell.graph.info.part.8 tuvx_no2.json tuvx_upper_atm.csv; do
        cp "$baseline/$name" "$case_dir/$name"
    done
    cp "$run_root/inert.yaml" "$case_dir/lnox_o3.yaml"
    ln -s "$init" "$case_dir/supercell_init.nc"
    ln -s "$exe" "$case_dir/atmosphere_model"
    ln -s "$baseline/data" "$case_dir/data"
    shopt -s nullglob
    for path in "$baseline"/*.TBL "$baseline"/*DATA*; do
        ln -s "$path" "$case_dir/$(basename "$path")"
    done
    shopt -u nullglob
    python - "$case_dir" "$duration" "${@:5}" <<'PY'
import re
import sys
from pathlib import Path
run, duration, options = Path(sys.argv[1]), sys.argv[2], dict(a.split('=', 1) for a in sys.argv[3:])
nml = run / 'namelist.atmosphere'
text = nml.read_text()

def setopt(name, value, group):
    global text
    pattern = rf"^(\s*){name}\s*=.*$"
    if re.search(pattern, text, flags=re.M):
        text = re.sub(pattern, rf"\g<1>{name} = {value}", text, count=1, flags=re.M)
    else:
        text, count = re.subn(rf"^&{group}\s*$", f"&{group}\n    {name} = {value}", text, count=1, flags=re.M)
        assert count == 1, f'no &{group} group'

setopt('config_run_duration', f"'{duration}'", 'nhyd_model')
setopt('config_pbl_scheme', f"'{options.get('pbl', 'bl_ysu')}'", 'physics')
setopt('config_sfclayer_scheme', f"'{options.get('sfclayer', 'sf_monin_obukhov')}'", 'physics')
if options.get('transport', 'default') != 'default':
    setopt('config_micm_pbl_transport', options['transport'], 'physics')
if 'pblmix' in options:
    setopt('config_ysu_pblmix', options['pblmix'], 'physics')
if 'pbl_interval' in options:
    setopt('config_pbl_interval', f"'{options['pbl_interval']}'", 'physics')
setopt('config_lnox_source_rate', '0.0', 'lnox')
setopt('config_tuvx_config_file', "''", 'photolysis')
setopt('config_tuvx_top_extension', '.false.', 'photolysis')
setopt('config_tuvx_upper_column_mode', "'none'", 'photolysis')
setopt('config_tuvx_extension_file', "''", 'photolysis')
setopt('config_j_no2_max', '0.0', 'photolysis')
nml.write_text(text)

streams = run / 'streams.atmosphere'
s = streams.read_text()
assert 'output_interval="00:10:00"' in s
streams.write_text(s.replace('output_interval="00:10:00"', 'output_interval="00:01:00"'))
listing = run / 'stream_list.atmosphere.output'
names = [n for n in listing.read_text().split() if n != 'j_jNO2']
names += ['hpbl', 'kpbl', 'hfx', 'ust', 'kzh']
listing.write_text('\n'.join(names) + '\n')
PY
}

status=0
fail() { echo "FAIL $*"; status=1; }

run_model() {
    # run_model CASE: returns the model's exit status
    local case_dir=$run_root/$1
    set +e
    (cd "$case_dir" && mpiexec -n 8 ./atmosphere_model > run.out 2>&1)
    local rc=$?
    set -e
    return $rc
}

has_log() {
    # has_log CASE TEXT: TEXT appears in the case's model logs. grep reads the files itself: a
    # pipe into grep -q fails under pipefail when grep exits before the writer finishes.
    local case_dir=$run_root/$1
    grep -qsF -- "$2" "$case_dir"/log.atmosphere.*.out "$case_dir"/log.atmosphere.*.err "$case_dir/run.out"
}

expect_stop() {
    # expect_stop CASE MESSAGE
    if run_model "$1"; then
        fail "$1: the model ran; expected it to stop"
    elif has_log "$1" "$2"; then
        echo "PASS $1: the model stopped with '$2'"
    else
        fail "$1: the model stopped without '$2'; see $run_root/$1"
    fi
}

echo "== gate cases"
stage gate-mynn "$run_root/stable_init.nc" "$model" 00:00:06 pbl=bl_mynn sfclayer=sf_mynn transport=.true.
expect_stop gate-mynn "MYNN gas transport is not yet qualified"
stage gate-interval "$run_root/stable_init.nc" "$model" 00:00:06 transport=.true. pbl_interval=00:00:06
expect_stop gate-interval "config_pbl_interval must be 'none'"
stage gate-no-sfclayer "$run_root/stable_init.nc" "$model" 00:00:06 sfclayer=off transport=.true.
expect_stop gate-no-sfclayer "wrong choice for surface layer scheme with YSU PBL"
stage gate-pbl-off "$run_root/stable_init.nc" "$model" 00:00:06 pbl=off sfclayer=off transport=.true.
if ! run_model gate-pbl-off; then
    fail "gate-pbl-off: the model stopped; see $run_root/gate-pbl-off"
elif has_log gate-pbl-off "PBL: resolved scheme = off" && \
     has_log gate-pbl-off "PBL gas transport enabled = F"; then
    echo "PASS gate-pbl-off: runs; PBL gas transport logged as not enabled"
else
    fail "gate-pbl-off: missing transport log lines"
fi

echo "== YSU runs (300 s, output every 60 s)"
stage stable-off "$run_root/stable_init.nc" "$model" 00:05:00
stage stable-on "$run_root/stable_init.nc" "$model" 00:05:00 transport=.true.
stage convective-off "$run_root/convective_init.nc" "$model" 00:05:00 transport=.false.
stage convective-on "$run_root/convective_init.nc" "$model" 00:05:00 transport=.true.
stage convective-topdown-off "$run_root/convective_init.nc" "$model" 00:05:00 transport=.false. pblmix=.true.
stage convective-topdown-on "$run_root/convective_init.nc" "$model" 00:05:00 transport=.true. pblmix=.true.
cases=(stable-off stable-on convective-off convective-on convective-topdown-off convective-topdown-on)
if [[ -n "$pre_p3_model" ]]; then
    stage stable-off-pre-p3 "$run_root/stable_init.nc" "$pre_p3_model" 00:05:00
    cases+=(stable-off-pre-p3)
fi
ran=1
for c in "${cases[@]}"; do
    run_model "$c" || { fail "$c: the model exited with an error; see $run_root/$c"; ran=0; }
done
for c in stable-on convective-on convective-topdown-on; do
    has_log "$c" "PBL gas transport enabled = T" || fail "$c: PBL gas transport not logged as enabled"
done
has_log stable-off "PBL gas transport enabled = F" || fail "stable-off: transport logged as enabled"

if (( ran )); then
    echo "== stable surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/stable-on" "$run_root/stable-off" \
        --json "$run_root/stable.json" || status=1
    echo "== convective surface layer"
    python "$here/check_pbl_ysu.py" mixing "$run_root/convective-on" "$run_root/convective-off" \
        --json "$run_root/convective.json" || status=1
    # NO2's step (about 1.2 km) lies above both PBL tops: both regimes mix it with YSU's
    # free-atmosphere coefficient, so only the in-PBL gradients are compared.
    python "$here/check_pbl_ysu.py" stronger "$run_root/convective.json" "$run_root/stable.json" \
        --gases qNO qO3 || status=1
    echo "== convective, top-down mixing on (config_ysu_pblmix)"
    python "$here/check_pbl_ysu.py" mixing "$run_root/convective-topdown-on" "$run_root/convective-topdown-off" \
        --json "$run_root/convective-topdown.json" || status=1
    echo "== transport off is unchanged from the pre-P3 executable"
    if [[ -n "$pre_p3_model" ]]; then
        python "$here/check_pbl_ysu.py" compare "$run_root/stable-off" "$run_root/stable-off-pre-p3" || status=1
    else
        echo "SKIP no PRE_P3_MODEL given"
    fi
fi

if (( status == 0 )); then
    echo "PASS: YSU gas PBL qualification"
else
    echo "FAIL: YSU gas PBL qualification"
fi
exit "$status"

# Shared helpers for the in-model PBL gas transport qualifications (check_pbl_ysu.sh,
# check_pbl_mynn.sh; docs/plans/pbl-transport.md P3 and P4). Source after setting:
#   baseline  staged supercell baseline directory (with environment.sh)
#   run_root  new, empty run directory
#   here      testing/chemistry_transport
#   PBL_DEFAULT, SFCLAYER_DEFAULT   scheme names used when a case does not set pbl= / sfclayer=
# Cases run the supercell with an inert mechanism (titration A = 0, j = 0, TUV-x off) and lightning
# off, so the gases are passive; gas profiles with vertical structure come from make_pbl_init.py.

status=0
fail() { echo "FAIL $*"; status=1; }

prepare_inputs() {
    # stable_init.nc (baseline surface, 288 K), convective_init.nc (skintemp = sst = 305 K), inert.yaml
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
}

stage() {
    # stage CASE INIT MODEL DURATION KEY=VALUE...  namelist options: pbl, sfclayer, transport
    # (default leaves config_micm_pbl_transport unset), pblmix, pbl_interval, mynn_edmf, mynn_edmf_dd,
    # mynn_edmf_output; extra_output adds comma-separated fields to the output stream.
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
    PBL_DEFAULT=$PBL_DEFAULT SFCLAYER_DEFAULT=$SFCLAYER_DEFAULT \
    python - "$case_dir" "$duration" "${@:5}" <<'PY'
import os
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
setopt('config_pbl_scheme', f"'{options.get('pbl', os.environ['PBL_DEFAULT'])}'", 'physics')
setopt('config_sfclayer_scheme', f"'{options.get('sfclayer', os.environ['SFCLAYER_DEFAULT'])}'", 'physics')
if options.get('transport', 'default') != 'default':
    setopt('config_micm_pbl_transport', options['transport'], 'physics')
if 'pblmix' in options:
    setopt('config_ysu_pblmix', options['pblmix'], 'physics')
if 'pbl_interval' in options:
    setopt('config_pbl_interval', f"'{options['pbl_interval']}'", 'physics')
for key in ('mynn_edmf', 'mynn_edmf_dd', 'mynn_edmf_output'):
    if key in options:
        setopt(f'config_{key}', options[key], 'physics')
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
names += [n for n in options.get('extra_output', '').split(',') if n]
listing.write_text('\n'.join(names) + '\n')
PY
}

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

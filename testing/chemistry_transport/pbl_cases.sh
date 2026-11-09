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
    # mynn_edmf_output, budget (config_micm_transport_budget), convection, chemistry_interval,
    # mechanism=PATH (a MICM file instead of the inert one), j_no2_max (default 0, as inert),
    # lnox_source_rate (default 0: lightning off); ranks=N runs
    # on N ranks (run_model) with supercell.graph.info.part.N from $CHEMPAS_PARTITIONS (default
    # ~/Data/MPAS/supercell); streams:
    # output_interval (history, default 00:01:00), restart_interval (restart stream; default as the
    # baseline); restart_from=PATH restarts from that restart file (config_do_restart, with the start time
    # from its name); extra_output adds comma-separated fields to the output stream.
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
if 'v_eddy_visc2' in options:
    setopt('config_v_mom_eddy_visc2', options['v_eddy_visc2'], 'nhyd_model')
    setopt('config_v_theta_eddy_visc2', options['v_eddy_visc2'], 'nhyd_model')
if 'budget' in options:
    setopt('config_micm_transport_budget', options['budget'], 'physics')
if 'convection' in options:
    setopt('config_convection_scheme', f"'{options['convection']}'", 'physics')
if 'chemistry_interval' in options:
    setopt('config_chemistry_interval', options['chemistry_interval'], 'chemistry')
if 'restart_from' in options:
    src = Path(options['restart_from']).resolve()
    stamp = src.name[len('restart.'):-len('.nc')]              # 0000-01-01_18.01.12
    day, clock = stamp.split('_')
    setopt('config_start_time', f"'{day}_{clock.replace('.', ':')}'", 'nhyd_model')
    setopt('config_do_restart', '.true.', 'restart')
    os.symlink(src, run / src.name)
setopt('config_lnox_source_rate', options.get('lnox_source_rate', '0.0'), 'lnox')
setopt('config_tuvx_config_file', "''", 'photolysis')
setopt('config_tuvx_top_extension', '.false.', 'photolysis')
setopt('config_tuvx_upper_column_mode', "'none'", 'photolysis')
setopt('config_tuvx_extension_file', "''", 'photolysis')
setopt('config_j_no2_max', options.get('j_no2_max', '0.0'), 'photolysis')
nml.write_text(text)
if 'ranks' in options:
    ranks = int(options['ranks'])
    (run / 'ranks').write_text(f'{ranks}\n')
    if ranks > 1:
        parts = Path(os.environ.get('CHEMPAS_PARTITIONS', Path.home() / 'Data/MPAS/supercell'))
        (run / f'supercell.graph.info.part.{ranks}').write_bytes((parts / f'supercell.graph.info.part.{ranks}').read_bytes())
if 'mechanism' in options:
    (run / 'lnox_o3.yaml').unlink()
    (run / 'lnox_o3.yaml').write_text(Path(options['mechanism']).read_text())

streams = run / 'streams.atmosphere'
s = streams.read_text()
assert 'output_interval="00:10:00"' in s
s = s.replace('output_interval="00:10:00"', f"output_interval=\"{options.get('output_interval', '00:01:00')}\"")
if 'restart_interval' in options:
    m = re.search(r'<immutable_stream\s+name="restart".*?>', s, flags=re.S)
    tag, n = re.subn(r'output_interval="[^"]*"', f'output_interval="{options["restart_interval"]}"', m.group(0))
    assert n == 1, 'restart stream has no output_interval'
    s = s[:m.start()] + tag + s[m.end():]
streams.write_text(s)
listing = run / 'stream_list.atmosphere.output'
names = [n for n in listing.read_text().split() if n != 'j_jNO2']
names += ['hpbl', 'kpbl', 'hfx', 'ust', 'kzh']
names += [n for n in options.get('extra_output', '').split(',') if n]
listing.write_text('\n'.join(names) + '\n')
PY
}

run_model() {
    # run_model CASE: returns the model's exit status. Eight ranks, or the count staged with ranks=N.
    local case_dir=$run_root/$1 ranks=8
    [[ -f "$case_dir/ranks" ]] && ranks=$(< "$case_dir/ranks")
    set +e
    (cd "$case_dir" && mpiexec -n "$ranks" ./atmosphere_model > run.out 2>&1)
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

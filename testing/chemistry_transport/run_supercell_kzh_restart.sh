#!/usr/bin/env bash
# Restart a finished supercell run (run_supercell_lnox_pbl.sh) from its restart
# file for a short time and write YSU's scalar diffusivity kzh, which the
# standard output stream does not carry. MPAS's kzh is exch_h on the bottom face
# of each layer; above the PBL top it is the coefficient the MICM gases see
# (bl_ysu.F90: xkzq = xkzh there). check_pbl_lnox_mixing.py --kzh-run reads it.
# The meteorology does not depend on gas transport, so either run of a pair
# gives the diffusivity of both.
#
# Usage: run_supercell_kzh_restart.sh RUN_DIR OUT_DIR [DURATION]
#   RUN_DIR   finished run with exactly one restart.*.nc (left unchanged)
#   OUT_DIR   new directory; output_kzh.nc (kzh, kpbl and the standard fields)
#             is written there once, after DURATION
#   DURATION  hh:mm:ss, default 00:00:30 (10 steps of 3 s)
set -euo pipefail
if (( $# < 2 || $# > 3 )); then
    echo "Usage: $0 RUN_DIR OUT_DIR [DURATION]" >&2
    exit 2
fi
run=$(cd "$1" && pwd)
out=$2
duration=${3:-00:00:30}
if [[ -e "$out" ]]; then
    echo "Use a new directory; existing runs are preserved: $out" >&2
    exit 2
fi
shopt -s nullglob
restarts=("$run"/restart.*.nc)
shopt -u nullglob
if (( ${#restarts[@]} != 1 )); then
    echo "Expected exactly one restart file in $run, found ${#restarts[@]}" >&2
    exit 2
fi
restart=${restarts[0]}
stamp=$(basename "$restart" .nc)
stamp=${stamp#restart.}                          # e.g. 0000-01-01_20.00.00
start="${stamp%%_*}_$(tr . : <<< "${stamp#*_}")"  # e.g. 0000-01-01_20:00:00

mkdir -p "$out/provenance"
out=$(cd "$out" && pwd)
for path in "$run"/*; do
    name=$(basename "$path")
    case $name in
        namelist.atmosphere|streams.atmosphere|stream_list.atmosphere.output|output*.nc|restart.*|\
        restart_timestamp|log.*|run.out|provenance|bin|atmosphere_model) ;;
        *) ln -s "$(readlink -f "$path")" "$out/$name" ;;
    esac
done
ln -s "$restart" "$out/$(basename "$restart")"
ln -s "$(readlink -f "$run/atmosphere_model")" "$out/atmosphere_model"

python3 - "$run" "$out" "$start" "$duration" <<'PY'
import re
import sys
from pathlib import Path
run, out, start, duration = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]

text = (run / "namelist.atmosphere").read_text()
for name, value in (("config_start_time", f"'{start}'"), ("config_run_duration", f"'{duration}'"),
                    ("config_do_restart", ".true.")):
    text, n = re.subn(rf"^(\s*){name}\s*=.*$", rf"\g<1>{name} = {value}", text, count=1, flags=re.M)
    assert n == 1, f"{name} not in namelist.atmosphere"
(out / "namelist.atmosphere").write_text(text)

streams = (run / "streams.atmosphere").read_text()
def edit_tag(s, stream, **attrs):
    """Set attributes in the opening tag of the named stream."""
    m = re.search(rf'<(?:immutable_)?stream\s+name="{stream}".*?>', s, flags=re.S)
    assert m, f"no {stream} stream"
    tag = m.group(0)
    for key, value in attrs.items():
        tag, n = re.subn(rf'{key}="[^"]*"', f'{key}="{value}"', tag)
        assert n == 1, f"{stream} stream has no {key}"
    return s[:m.start()] + tag + s[m.end():]
streams = edit_tag(streams, "restart", output_interval="none")
streams = edit_tag(streams, "output", filename_template="output_kzh.nc", output_interval=duration)
(out / "streams.atmosphere").write_text(streams)

names = (run / "stream_list.atmosphere.output").read_text().split()
names += [n for n in ("hpbl", "kzh", "kpbl") if n not in names]
(out / "stream_list.atmosphere.output").write_text("\n".join(names) + "\n")
PY

source "$out/environment.sh"
printf 'restart of %s at %s for %s, writing kzh\n' "$run" "$start" "$duration" > "$out/provenance/purpose"
sha256sum "$(readlink -f "$out/atmosphere_model")" > "$out/provenance/model-sha256.txt"
cd "$out"
set +e
mpiexec -n 8 ./atmosphere_model > run.out 2>&1
status=$?
set -e
echo "$status" > provenance/run-exit-status
exit "$status"

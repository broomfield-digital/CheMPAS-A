#!/usr/bin/env bash
# Run the standard supercell LNOx case (the 2-h baseline configuration: initial
# NO = NO2 = 0, O3 = 50 ppbv, altitude-gated lightning NOx, 10-min output) with
# a PBL scheme, and MICM gas PBL transport off or on. A pair of runs that differ
# only in TRANSPORT isolates the effect of gas PBL transport
# (docs/plans/pbl-transport.md); the baseline itself has PBL off.
#
# Usage: run_supercell_lnox_pbl.sh BASELINE_DIR RUN_DIR off|on [MODEL [PBL]]
#   PBL: ysu (default; with sf_monin_obukhov) or mynn (sf_mynn, MPAS defaults)
set -euo pipefail
if (( $# < 3 || $# > 5 )); then
    echo "Usage: $0 BASELINE_DIR RUN_DIR off|on [MODEL [ysu|mynn]]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
run_dir=$2
case "$3" in on) transport=.true. ;; off) transport=.false. ;; *) echo "TRANSPORT must be on or off" >&2; exit 2 ;; esac
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${4:-$repo/atmosphere_model}")
case "${5:-ysu}" in
    ysu) pbl=bl_ysu; sfclayer=sf_monin_obukhov ;;
    mynn) pbl=bl_mynn; sfclayer=sf_mynn ;;
    *) echo "PBL must be ysu or mynn" >&2; exit 2 ;;
esac
test -x "$model"
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
shopt -u nullglob

python3 - "$run_dir" "$pbl" "$sfclayer" "$transport" <<'PY'
import re
import sys
from pathlib import Path
run, pbl, sfclayer, transport = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
nml = run / "namelist.atmosphere"
text = nml.read_text()

def setopt(name, value, group):
    global text
    pattern = rf"^(\s*){name}\s*=.*$"
    if re.search(pattern, text, flags=re.M):
        text = re.sub(pattern, rf"\g<1>{name} = {value}", text, count=1, flags=re.M)
    else:
        text, count = re.subn(rf"^&{group}\s*$", f"&{group}\n    {name} = {value}", text, count=1, flags=re.M)
        assert count == 1, f"no &{group} group"

setopt("config_pbl_scheme", f"'{pbl}'", "physics")
setopt("config_sfclayer_scheme", f"'{sfclayer}'", "physics")
setopt("config_micm_pbl_transport", transport, "physics")
nml.write_text(text)
listing = run / "stream_list.atmosphere.output"
names = listing.read_text().split()
names += [n for n in ("hpbl", "hfx") if n not in names]
listing.write_text("\n".join(names) + "\n")
PY

source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_dir/provenance/chempas_revision"
sha256sum "$model" > "$run_dir/provenance/model-sha256.txt"
printf 'pbl=%s sfclayer=%s transport=%s\n' "$pbl" "$sfclayer" "$transport" > "$run_dir/provenance/case"
cd "$run_dir"
date -u +%FT%TZ > provenance/run-started-utc
set +e
mpiexec -n 8 ./atmosphere_model > run.out 2>&1
status=$?
set -e
echo "$status" > provenance/run-exit-status
date -u +%FT%TZ > provenance/run-finished-utc
exit "$status"

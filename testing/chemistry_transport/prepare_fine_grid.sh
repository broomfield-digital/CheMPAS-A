#!/usr/bin/env bash
# Prepare the finer-vertical-grid supercell fixture (P6 of
# docs/plans/pbl-transport.md: "a finer-grid fixture"). The baseline grid has
# interfaces at 0, 299, 712 and 1182 m, so its PBL spans one or two layers;
# this one has 50-m layers to 1.5 km (make_fine_zeta_levels.py) under the same
# 50-km top, on the same horizontal mesh and partition.
#
# Usage: prepare_fine_grid.sh BASELINE_DIR OUT_DIR
#
# Runs the baseline's init_atmosphere_model (supercell, config_init_case = 5)
# on eight ranks with the new levels, adds the LNOx gases (init_lnox_o3.py),
# then writes the PBL gas profiles of make_pbl_init.py:
#   stable_init.nc       baseline surface (288 K)
#   convective_init.nc   skintemp = sst = 305 K
# OUT_DIR must not exist. check_pbl_fine.sh runs the cases. Model runs on this
# grid need the supercell's explicit vertical eddy viscosity off at dt = 3 s
# (500 m2/s on 50-m layers is unstable); check_pbl_fine.sh sets it to zero.
set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: $0 BASELINE_DIR OUT_DIR" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
out=$2
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
if [[ -e "$out" ]]; then
    echo "Use a new directory; existing fixtures are preserved: $out" >&2
    exit 2
fi
mkdir -p "$out/provenance"
out=$(cd "$out" && pwd)
source "$baseline/environment.sh"

python "$here/make_fine_zeta_levels.py" "$out/supercell_zeta_levels.txt"
nlev=$(( $(wc -l < "$out/supercell_zeta_levels.txt") - 1 ))
for name in namelist.init_atmosphere streams.init_atmosphere supercell.graph.info.part.8; do
    cp "$baseline/$name" "$out/$name"
done
ln -s "$baseline/supercell_grid.nc" "$out/supercell_grid.nc"
cp "$baseline/init_atmosphere_model" "$out/init_atmosphere_model"
python3 - "$out/namelist.init_atmosphere" "$nlev" <<'PY'
import re
import sys
from pathlib import Path
path, nlev = Path(sys.argv[1]), sys.argv[2]
text, n = re.subn(r"^(\s*config_nvertlevels\s*=).*$", rf"\g<1> {nlev}", path.read_text(), count=1, flags=re.M)
assert n == 1, "config_nvertlevels not in namelist.init_atmosphere"
assert "config_specified_zeta_levels = 'supercell_zeta_levels.txt'" in text
path.write_text(text)
PY

cd "$out"
mpiexec -n 8 ./init_atmosphere_model > init_run.out 2>&1
test -s supercell_init.nc
python "$repo/scripts/init_lnox_o3.py" -i supercell_init.nc > initialize_tracers.log 2>&1
python "$here/make_pbl_init.py" supercell_init.nc stable_init.nc
python "$here/make_pbl_init.py" supercell_init.nc convective_init.nc --skin-temperature 305
git -C "$repo" rev-parse HEAD > provenance/chempas_revision
sha256sum supercell_zeta_levels.txt supercell_init.nc stable_init.nc convective_init.nc \
    "$baseline/supercell_grid.nc" init_atmosphere_model > provenance/sha256.txt
echo "fine grid: $nlev layers in $out"

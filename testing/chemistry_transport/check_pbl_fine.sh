#!/usr/bin/env bash
# In-model PBL gas mixing on the finer vertical grid (P6 of
# docs/plans/pbl-transport.md). On the baseline grid the PBL spans one or two
# layers and MYNN's mass-flux plumes never start; on this grid (50-m layers to
# 1.5 km, prepare_fine_grid.sh) the convective PBL spans many layers.
#
# Usage: check_pbl_fine.sh BASELINE_DIR FINE_DIR RUN_ROOT [MODEL]
#
# Cases, 10 min each, convective init (skintemp = sst = 305 K), inert mechanism,
# transport off and on, budget logging on with transport, and the supercell's
# explicit vertical eddy viscosity off (config_v_mom_eddy_visc2 =
# config_v_theta_eddy_visc2 = 0): the PBL scheme mixes the resolved boundary
# layer, and the explicit term, applied once per step, needs K*dt/dz**2 <= 0.5,
# which 500 m2/s on 50-m layers at dt = 3 s (0.6) breaks within two minutes:
#   ysu-off/on               YSU
#   mynn-edmf-off/on         MYNN with its updraft plumes (config_mynn_edmf = 1),
#                            edmf_a written (config_mynn_edmf_output = 1)
#   mynn-9-off/on            the same at a 9-s PBL interval
# Checks: check_pbl_ysu.py mixing for each pair; check_pbl_fine.py mixed (the
# PBL spans at least 8 layers; NO inside it mixes at the rate of the scheme's own
# diffusivity, within a factor of 2 of an offline replay with the output kzh; for
# MYNN the plumes are active); check_transport_budget.py on each run with
# transport. Ten minutes is shorter than the PBL mixing time (about 16 min for
# YSU here, longer for MYNN's cold-started TKE), so the PBL is not yet well mixed.
# RUN_ROOT must not exist.
set -euo pipefail
if (( $# < 3 || $# > 4 )); then
    echo "Usage: $0 BASELINE_DIR FINE_DIR RUN_ROOT [MODEL]" >&2
    exit 2
fi
baseline=$(cd "$1" && pwd)
fine=$(cd "$2" && pwd)
run_root=$3
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
model=$(realpath "${4:-$repo/atmosphere_model}")
test -x "$model"
test -s "$fine/convective_init.nc"
if [[ -e "$run_root" ]]; then
    echo "Use a new run root; existing runs are preserved: $run_root" >&2
    exit 2
fi
mkdir -p "$run_root/provenance"
run_root=$(cd "$run_root" && pwd)
source "$baseline/environment.sh"
git -C "$repo" rev-parse HEAD > "$run_root/provenance/chempas_revision"
git -C "$repo" status --short > "$run_root/provenance/source-status"
sha256sum "$model" "$fine/convective_init.nc" > "$run_root/provenance/sha256.txt"

PBL_DEFAULT=bl_ysu
SFCLAYER_DEFAULT=sf_monin_obukhov
# shellcheck source=pbl_cases.sh
source "$here/pbl_cases.sh"
prepare_inputs    # the inert mechanism (and baseline-grid inits, unused here)

init=$fine/convective_init.nc
dynamics=(v_eddy_visc2=0.0)
mynn=(pbl=bl_mynn sfclayer=sf_mynn mynn_edmf=1 mynn_edmf_output=1 extra_output=edmf_a)
stage ysu-off "$init" "$model" 00:10:00 "${dynamics[@]}" transport=.false.
stage ysu-on "$init" "$model" 00:10:00 "${dynamics[@]}" transport=.true. budget=.true.
stage mynn-edmf-off "$init" "$model" 00:10:00 "${dynamics[@]}" "${mynn[@]}" transport=.false.
stage mynn-edmf-on "$init" "$model" 00:10:00 "${dynamics[@]}" "${mynn[@]}" transport=.true. budget=.true.
stage mynn-9-off "$init" "$model" 00:10:00 "${dynamics[@]}" "${mynn[@]}" pbl_interval=00:00:09 transport=.false.
stage mynn-9-on "$init" "$model" 00:10:00 "${dynamics[@]}" "${mynn[@]}" pbl_interval=00:00:09 transport=.true. budget=.true.
cases=(ysu-off ysu-on mynn-edmf-off mynn-edmf-on mynn-9-off mynn-9-on)
ran=1
for (( n = 0; n < ${#cases[@]}; n += 2 )); do
    pids=()
    for c in "${cases[@]:n:2}"; do
        run_model "$c" & pids+=("$!")
    done
    m=0
    for c in "${cases[@]:n:2}"; do
        wait "${pids[m]}" || { fail "$c: the model exited with an error; see $run_root/$c"; ran=0; }
        m=$((m + 1))
    done
done

if (( ran )); then
    for pair in ysu:3: mynn-edmf:3:--plumes mynn-9:9:--plumes; do
        IFS=: read -r name interval plumes <<< "$pair"
        echo "== $name"
        python "$here/check_pbl_ysu.py" mixing "$run_root/$name-on" "$run_root/$name-off" \
            --json "$run_root/$name-mixing.json" || status=1
        # shellcheck disable=SC2086
        python "$here/check_pbl_fine.py" mixed "$run_root/$name-on" "$run_root/$name-off" $plumes \
            --json "$run_root/$name-fine.json" || status=1
        python "$here/check_transport_budget.py" "$run_root/$name-on/log.atmosphere.0000.out" --expect-gases 3 \
            --expect-steps 200 --dt 3 --interval "$interval" --json "$run_root/$name-budget.json" || status=1
    done
fi

if (( status == 0 )); then
    echo "PASS: PBL gas mixing on the finer grid"
else
    echo "FAIL: PBL gas mixing on the finer grid"
fi
exit "$status"

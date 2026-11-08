#!/usr/bin/env bash
# Figures for P2 of docs/plans/pbl-transport.md: the column diffusion kernel and
# the realized-increment proof.
#
# Usage: make_p2_plots.sh REALIZED_RUN_ROOT OUTDIR
#
# REALIZED_RUN_ROOT is a RUN_ROOT written by ../check_realized_increment.sh.
# Builds dump_pbl_diffusion.F90 against the production kernel with the flags in
# .build_opts.atmosphere (build the atmosphere core first, and run in the same
# environment), writes the kernel data to a temporary directory, then draws
# p2_1..p2_4 into OUTDIR with plot_p2_evidence.py.
set -euo pipefail
if (( $# != 2 )); then
    echo "Usage: $0 REALIZED_RUN_ROOT OUTDIR" >&2
    exit 2
fi
run_root=$(cd "$1" && pwd)
mkdir -p "$2"
outdir=$(cd "$2" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
opts=$repo/.build_opts.atmosphere
test -f "$opts" || { echo "Missing $opts; build the atmosphere core first." >&2; exit 2; }
opt() { sed -n "s/^$1=//p" "$opts"; }
FC=$(opt FC)
FFLAGS=$(opt FFLAGS)
CPPFLAGS=$(opt CPPFLAGS)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir "$work/data"
(
    cd "$work"
    # shellcheck disable=SC2086
    $FC $CPPFLAGS $FFLAGS -c "$repo/src/core_atmosphere/physics/mpas_atmphys_pbl_diffusion.F" -I"$repo/src/framework"
    # shellcheck disable=SC2086
    $FC $CPPFLAGS $FFLAGS -c "$here/dump_pbl_diffusion.F90" -I"$repo/src/framework" -I.
    # shellcheck disable=SC2086
    $FC $FFLAGS -o dump_pbl_diffusion dump_pbl_diffusion.o mpas_atmphys_pbl_diffusion.o
    ./dump_pbl_diffusion "$repo/testing/chemistry_transport/columns" "$work/data"
)
python "$here/plot_p2_evidence.py" "$work/data" "$run_root" "$outdir"

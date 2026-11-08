#!/usr/bin/env bash
# Build and run the standalone Fortran unit tests for the MICM gas transport
# foundation (docs/plans/pbl-transport.md, P1).
#
# Usage: run_unit_tests.sh [TEST_NAME ...]    (default: every test_*.F90 here)
#
# Compiler, flags and libraries come from the repository's
# .build_opts.atmosphere, so build the atmosphere core once first and run this
# in the same environment (for example, source the baseline environment.sh).
# Production sources are compiled from src/ into a temporary directory; the
# build tree is not modified.
#
# A test may replace a production module with a test double: a file
# doubles/<TEST_NAME>/<production basename>.F90 is compiled instead of the
# production source of the same basename, for that test only.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
opts=$repo/.build_opts.atmosphere
if [[ ! -f "$opts" ]]; then
    echo "Missing $opts; build the atmosphere core first." >&2
    exit 2
fi
opt() { sed -n "s/^$1=//p" "$opts"; }
FC=$(opt FC)
FFLAGS=$(opt FFLAGS)
CPPFLAGS=$(opt CPPFLAGS)
CPPINCLUDES=$(opt CPPINCLUDES)
LIBS=$(opt LIBS)
command -v "$FC" > /dev/null || { echo "$FC is not on PATH" >&2; exit 2; }
# With gfortran, unit tests also run with runtime checks (bounds, unallocated
# arrays, pointers); later options override the production -O3.
if "$FC" --version 2>/dev/null | grep -q "GNU Fortran"; then
    FFLAGS="$FFLAGS -O0 -g -fcheck=all -fbacktrace"
fi

framework=$repo/src/framework
core=$repo/src/core_atmosphere
# The YSU and MYNN column tests need the MMM-physics export patches (idempotent).
"$core/physics/apply_mmm_patches.sh" > /dev/null
esmf_time=$repo/src/external/esmf_time_f90

# Production sources under test, in dependency order.
sources=(
    "$core/chemistry/musica/mpas_chem_transport_map.F"
    "$core/physics/mpas_atmphys_utilities.F"
    "$core/physics/ccpp_kind_types.F"
    "$core/physics/physics_mmm/bl_ysu.F90"
    "$core/physics/physics_wrf/module_bl_ysu.F"
    "$core/physics/physics_mmm/mynn_shared.F90"
    "$core/physics/physics_mmm/bl_mynn_subroutines.F90"
    "$core/physics/physics_mmm/bl_mynn.F90"
    "$core/physics/physics_wrf/bl_mynn_pre.F"
    "$core/physics/physics_wrf/bl_mynn_post.F"
    "$core/physics/physics_wrf/module_bl_mynn.F"
    "$core/physics/mpas_atmphys_tracer_column.F"
    "$core/physics/mpas_atmphys_pbl_diffusion.F"
    "$core/physics/mpas_atmphys_pbl_edmf.F"
    "$core/physics/mpas_atmphys_pbl_gas.F"
    "$core/physics/mpas_atmphys_cu_gas.F"
    "$core/physics/mpas_atmphys_tracer_transport.F"
)
support=("$here/unit_assert.F90")
# Test helpers that use production modules, compiled after the sources.
late_support=("$here/mynn_column.F90")
# Reference column fixtures (columns/make_diffusion_fixtures.py).
export CHEMPAS_UNIT_COLUMNS=${CHEMPAS_UNIT_COLUMNS:-$here/../columns}

if (( $# > 0 )); then
    tests=("$@")
else
    tests=()
    for path in "$here"/test_*.F90; do
        tests+=("$(basename "$path" .F90)")
    done
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

compile() {
    # shellcheck disable=SC2086
    $FC $CPPFLAGS $FFLAGS -c "$1" $CPPINCLUDES -I"$framework" -I.
}

failed=0
for name in "${tests[@]}"; do
    mkdir "$work/$name"
    cd "$work/$name"
    for src in "${support[@]}"; do
        compile "$src"
    done
    for src in "${sources[@]}"; do
        double=$here/doubles/$name/$(basename "$src" .F).F90
        if [[ -f "$double" ]]; then
            compile "$double"
        else
            compile "$src"
        fi
    done
    for src in "${late_support[@]}"; do
        compile "$src"
    done
    objects=(*.o)
    compile "$here/$name.F90"
    # shellcheck disable=SC2086
    $FC $FFLAGS -o "$name" "$name.o" "${objects[@]}" \
        "$framework/libframework.a" "$esmf_time/libesmf_time.a" $LIBS
    if ./"$name"; then
        echo "== $name: PASS"
    else
        echo "== $name: FAIL"
        failed=1
    fi
done
exit "$failed"

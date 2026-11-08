#!/usr/bin/env bash
# Self-test for src/core_atmosphere/physics/apply_mmm_patches.sh: a pristine
# MMM-physics file is patched to the recorded result, a patched file is left
# alone, and a file with other edits is refused. --revert (run by the physics
# Makefile before checkout_externals, which skips every external while one has
# local changes) restores the pinned file and a clean external; applying again
# restores the patched file with its earlier modification time, so an unchanged
# external does not trigger recompilation. Runs on copies in a temporary
# directory (the copy of physics_mmm is a git checkout of the pinned file); the
# build tree is not touched.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
physics=$repo/src/core_atmosphere/physics
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

failed=0
check() {
    if eval "$1"; then echo "PASS $2"; else echo "FAIL $2"; failed=1; fi
}

mkdir -p "$work/physics/physics_mmm"
cp "$physics/apply_mmm_patches.sh" "$work/physics/"
cp -r "$physics/patches" "$work/physics/"
mmm=$work/physics/physics_mmm
git -C "$physics/physics_mmm" show HEAD:bl_ysu.F90 > "$mmm/bl_ysu.F90"
git -C "$mmm" init -q
git -C "$mmm" add bl_ysu.F90
git -C "$mmm" -c user.name=test -c user.email=test@example.invalid commit -q -m pinned
manifest=$work/physics/patches/mmm-physics/MANIFEST
pristine=$(awk '$1 !~ /^#/ && $2 == "physics_mmm/bl_ysu.F90" {print $3}' "$manifest")
patched=$(awk '$1 !~ /^#/ && $2 == "physics_mmm/bl_ysu.F90" {print $4}' "$manifest")
script=$work/physics/apply_mmm_patches.sh
hash() { sha256sum "$mmm/bl_ysu.F90" | cut -d" " -f1; }

check '"$work/physics/apply_mmm_patches.sh" > "$work/out1" 2>&1' "pristine file is patched"
check '[[ $(sha256sum "$work/physics/physics_mmm/bl_ysu.F90" | cut -d" " -f1) == "$patched" ]]' \
      "patched file has the recorded hash"
check 'grep -q xkzqx "$work/physics/physics_mmm/bl_ysu.F90"' "export argument is present"
check '"$work/physics/apply_mmm_patches.sh" > "$work/out2" 2>&1 && grep -q "already patched" "$work/out2"' \
      "second application is a no-op"

touch -d "2001-02-03 04:05:06" "$mmm/bl_ysu.F90"
before=$(stat -c %Y "$mmm/bl_ysu.F90")
check '"$script" --revert > "$work/out4" 2>&1' "revert succeeds"
check '[[ $(hash) == "$pristine" ]]' "revert restores the pinned file"
check '[[ -z $(git -C "$mmm" status --porcelain) ]]' "revert leaves the external clean for checkout_externals"
check '"$script" --revert > "$work/out5" 2>&1 && grep -q "not patched" "$work/out5" && [[ $(hash) == "$pristine" ]]' \
      "reverting the pinned file is a no-op"
check '"$script" > "$work/out6" 2>&1 && [[ $(hash) == "$patched" ]]' "applying after a revert patches again"
check '[[ $(stat -c %Y "$mmm/bl_ysu.F90") == "$before" ]]' \
      "re-applied file keeps its earlier modification time"
check '[[ -z $(git -C "$mmm" status --porcelain --ignored | grep -v bl_ysu.F90) ]]' \
      "no state files are left in the external"

echo "! unrelated local edit" >> "$work/physics/physics_mmm/bl_ysu.F90"
check '! "$script" --revert > "$work/out7" 2>&1 && grep -q "unexpected contents" "$work/out7"' \
      "revert refuses a file with other edits"
check '! "$work/physics/apply_mmm_patches.sh" > "$work/out3" 2>&1 && grep -q "unexpected contents" "$work/out3"' \
      "a file with other edits is refused"

if (( failed )); then echo "test_apply_mmm_patches: FAIL"; exit 1; fi
echo "test_apply_mmm_patches: PASS"

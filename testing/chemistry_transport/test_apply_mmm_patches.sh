#!/usr/bin/env bash
# Self-test for src/core_atmosphere/physics/apply_mmm_patches.sh, over every
# MANIFEST entry (YSU, P3; MYNN, P4): a pristine MMM-physics file is patched
# to the recorded result, a patched file is left
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
manifest=$work/physics/patches/mmm-physics/MANIFEST
# every target file in MANIFEST, from the pinned external
mapfile -t targets < <(awk '$1 !~ /^#/ && NF == 4 {sub("^physics_mmm/", "", $2); print $2}' "$manifest")
for t in "${targets[@]}"; do
    git -C "$physics/physics_mmm" show "HEAD:$t" > "$mmm/$t"
done
git -C "$mmm" init -q
git -C "$mmm" add "${targets[@]}"
git -C "$mmm" -c user.name=test -c user.email=test@example.invalid commit -q -m pinned
script=$work/physics/apply_mmm_patches.sh
hash() { sha256sum "$mmm/$1" | cut -d" " -f1; }
recorded() { awk -v t="physics_mmm/$1" -v c="$2" '$1 !~ /^#/ && $2 == t {print $c}' "$manifest"; }
all_hashes() {
    # all_hashes COLUMN: every target has the MANIFEST hash in COLUMN (3 pinned, 4 patched)
    local t
    for t in "${targets[@]}"; do [[ $(hash "$t") == $(recorded "$t" "$1") ]] || return 1; done
}
check '(( ${#targets[@]} == 3 ))' "MANIFEST lists the YSU and the two MYNN targets"

check '"$script" > "$work/out1" 2>&1' "pristine files are patched"
check 'all_hashes 4' "patched files have the recorded hashes"
check 'grep -q xkzqx "$mmm/bl_ysu.F90"' "YSU export argument is present"
check 'grep -q "aw_up,ent_up" "$mmm/bl_mynn_subroutines.F90" && grep -q edmf_aw_up "$mmm/bl_mynn.F90"' \
      "MYNN plume export arguments are present"
check '"$script" > "$work/out2" 2>&1 && (( $(grep -c "already patched" "$work/out2") == 3 ))' \
      "second application is a no-op"

for t in "${targets[@]}"; do touch -d "2001-02-03 04:05:06" "$mmm/$t"; done
before=$(stat -c %Y "$mmm/bl_ysu.F90")
check '"$script" --revert > "$work/out4" 2>&1' "revert succeeds"
check 'all_hashes 3' "revert restores the pinned files"
check '[[ -z $(git -C "$mmm" status --porcelain) ]]' "revert leaves the external clean for checkout_externals"
check '"$script" --revert > "$work/out5" 2>&1 && (( $(grep -c "not patched" "$work/out5") == 3 )) && all_hashes 3' \
      "reverting the pinned files is a no-op"
check '"$script" > "$work/out6" 2>&1 && all_hashes 4' "applying after a revert patches again"
check 'for t in "${targets[@]}"; do [[ $(stat -c %Y "$mmm/$t") == "$before" ]] || exit 1; done' \
      "re-applied files keep their earlier modification times"
check '[[ -z $(git -C "$mmm" status --porcelain --ignored | grep -v "\.F90$") ]]' \
      "no state files are left in the external"

echo "! unrelated local edit" >> "$mmm/bl_mynn.F90"
check '! "$script" --revert > "$work/out7" 2>&1 && grep -q "unexpected contents" "$work/out7"' \
      "revert refuses a file with other edits"
check '! "$script" > "$work/out3" 2>&1 && grep -q "unexpected contents" "$work/out3"' \
      "a file with other edits is refused"

if (( failed )); then echo "test_apply_mmm_patches: FAIL"; exit 1; fi
echo "test_apply_mmm_patches: PASS"

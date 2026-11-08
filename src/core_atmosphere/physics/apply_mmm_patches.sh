#!/usr/bin/env bash
# Apply the CheMPAS-A patches in patches/mmm-physics/ to the MMM-physics external
# (see patches/mmm-physics/README.md). Each target must be the pinned original (it
# is patched and checked) or already patched (it is left alone); anything else stops
# the build.
#
# --revert restores the pinned originals. The physics Makefile runs it before
# checkout_externals, which skips every external while one has local changes, and
# applies the patches again afterwards. A reverted file is re-patched with its
# earlier modification time when the result is the same file, so an unchanged
# external is not recompiled; the time is kept in the external's git directory.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
patches=$here/patches/mmm-physics
mode=apply
if (( $# == 1 )) && [[ $1 == --revert ]]; then
    mode=revert
elif (( $# != 0 )); then
    echo "Usage: $0 [--revert]" >&2
    exit 2
fi

sha256() {
    if command -v sha256sum > /dev/null; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

mtime_record() {
    # mtime_record TARGET_FILE: reference file holding the patched hash and modification time
    local gitdir
    gitdir=$(git -C "$(dirname "$1")" rev-parse --absolute-git-dir 2> /dev/null) || return 1
    echo "$gitdir/chempas-mmm-patch.$(basename "$1")"
}

unexpected() {
    # unexpected TARGET SUM PATCH_FILE
    echo "apply_mmm_patches: $1 has unexpected contents (sha256 $2)." >&2
    echo "  Expected the pinned MMM-physics file or the result of $3." >&2
    echo "  To discard local edits: git -C $here/$(dirname "$1") checkout -- $(basename "$1")" >&2
    echo "  After an MMM-physics update, refresh the patch and patches/mmm-physics/MANIFEST." >&2
    exit 1
}

entries=()
while read -r patch_file target pristine patched || [[ -n ${patch_file:-} ]]; do
    [[ -z ${patch_file:-} || $patch_file == \#* ]] && continue
    entries+=("$patch_file $target $pristine $patched")
done < "$patches/MANIFEST"

if [[ $mode == revert ]]; then
    for (( n = ${#entries[@]} - 1; n >= 0; n-- )); do
        read -r patch_file target pristine patched <<< "${entries[n]}"
        file=$here/$target
        [[ -f $file ]] || continue
        sum=$(sha256 "$file")
        if [[ $sum == "$pristine" ]]; then
            echo "apply_mmm_patches: $target not patched"
            continue
        fi
        [[ $sum == "$patched" ]] || unexpected "$target" "$sum" "$patch_file"
        if record=$(mtime_record "$file"); then
            echo "$patched" > "$record"
            touch -r "$file" "$record"
        fi
        patch --quiet --reverse -p1 -d "$here" < "$patches/$patch_file"
        sum=$(sha256 "$file")
        if [[ $sum != "$pristine" ]]; then
            echo "apply_mmm_patches: reverting $patch_file did not restore the pinned $target" >&2
            exit 1
        fi
        echo "apply_mmm_patches: reverted $patch_file from $target"
    done
    exit 0
fi

for entry in "${entries[@]}"; do
    read -r patch_file target pristine patched <<< "$entry"
    file=$here/$target
    if [[ ! -f $file ]]; then
        echo "apply_mmm_patches: missing $target; run checkout_externals (the physics Makefile does)" >&2
        exit 1
    fi
    record=$(mtime_record "$file") || record=
    sum=$(sha256 "$file")
    if [[ $sum == "$patched" ]]; then
        echo "apply_mmm_patches: $target already patched ($patch_file)"
        [[ -z $record ]] || rm -f "$record"
        continue
    fi
    [[ $sum == "$pristine" ]] || unexpected "$target" "$sum" "$patch_file"
    patch --quiet --forward -p1 -d "$here" < "$patches/$patch_file"
    sum=$(sha256 "$file")
    if [[ $sum != "$patched" ]]; then
        patch --quiet --reverse -p1 -d "$here" < "$patches/$patch_file" || true
        echo "apply_mmm_patches: $target differs from the recorded result after $patch_file" >&2
        exit 1
    fi
    if [[ -n $record && -f $record ]]; then
        [[ $(cat "$record") == "$patched" ]] && touch -r "$record" "$file"
        rm -f "$record"
    fi
    echo "apply_mmm_patches: applied $patch_file to $target"
done

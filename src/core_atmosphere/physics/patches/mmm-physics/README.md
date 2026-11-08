# CheMPAS-A patches to MMM-physics

MMM-physics is an external pinned in `src/core_atmosphere/Externals.cfg`. The
pin is tag `20250616-MPASv8.3`, revision
`a4baf7f3243d1db0dbc5f63473f895bdbdc05c30`. `checkout_externals` fetches it
into `physics_mmm/`.

The MICM gas PBL transport (`docs/plans/pbl-transport.md`, P3) needs one
quantity that MMM-physics computes but does not export. The patches here add
that export. `../../apply_mmm_patches.sh` applies them, and the physics Makefile
runs it before it builds MMM-physics.

`checkout_externals` skips every external, with only a warning, while any
external has local changes. The Makefile therefore runs
`apply_mmm_patches.sh --revert` before `checkout_externals`, which restores the
pinned files and a clean external, and applies the patches again afterwards.
When an external has not changed, the re-patched file keeps its earlier
modification time, so it is not recompiled. The script keeps that time in the
external's git directory.

The script checks each target file's SHA-256 against `MANIFEST`:

- the pinned original: the patch is applied, and the result is checked
  against the recorded patched hash;
- the patched hash: nothing is done;
- anything else: it stops, and prints how to discard local edits.

An MMM-physics update in `Externals.cfg` is therefore checked out, and the
build then fails loudly until the patch and `MANIFEST` are refreshed.
The preferred long-term route is the same export upstream, in a versioned
MMM-physics release recorded in `Externals.cfg`. The patches can then be
removed.

| Patch | Change |
| --- | --- |
| `bl_ysu-export-xkzq.patch` | Adds the optional `intent(out)` argument `xkzqx(its:,:)` to `bl_ysu_run`. It holds the exchange coefficient [m² s⁻¹] that YSU's moisture and passive-tracer (`qmix`) solve uses at the interface between layers k and k+1: `xkzh` above the PBL top, and the entrainment-zone adjustment applied. When the argument is absent, results are unchanged. |

MMM-physics is distributed under the BSD 3-Clause License, Copyright (c) 2022,
National Center for Atmospheric Research. The patches contain portions of its
source. That license applies to them, and its text is in `LICENSE` here (copied
from the pinned revision). The patched files keep their original notices.

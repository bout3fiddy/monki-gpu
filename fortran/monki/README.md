# MONKI reference

MONKI 1.0.1 by Victor J.H. Trees, supplied as `September2026/MC_devel30` without
a Git revision. `upstream/` retains 22 byte-identical files: twelve Fortran sources,
three inputs, three archived NetCDF outputs, three optical tables and the MIT license.
These support the 340 nm clear, layered-cloud and box-cloud fixtures.

`manifest.json` records original SHA-256 hashes, the previous manifest hash and
inventory reduction. `scripts/build_reference.py` checks inventory and hashes;
generated dimensions and binaries go in `build/`. Manuals are in `docs/monki/`.

`PROGRAM/modules/ncdf_utilities_mod.f90` retains its LGPL-3.0-or-later notice;
MONKI's MIT license does not override it. All original notices remain intact.
Archived outputs are never regenerated. Their compiler/RNG details are unknown,
so comparisons are statistical; these fixtures do not validate other scenes.

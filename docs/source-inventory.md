# Source inventory

This document records the source inventory for `sharedLibs`.

It currently combines confirmed project information with an operational inventory plan. The repository is not yet fully indexed for code search, so this document should be refined as soon as the actual source tree is inspected locally or through a complete file tree listing.

## Confirmed repository components

The original project history identifies four main components:

| Component | Category | Origin | Current role | Development status | Recommended action |
| --- | --- | --- | --- | --- | --- |
| `libmisc` | Core | Author-developed / CPTEC-related | General-purpose utility library. | Active target. | Consolidate, document, build, and test first. |
| `libax` | Core | Author-developed / CPTEC-related | Scientific and meteorological data access library. | Active target. | Document API, dependencies, supported formats, and `w3lib` usage. |
| `sigioBAM` | Legacy | BAM-related local code | Spectral-file routines for BAM workflows. | Maintained legacy. | Review routines and migrate useful code when appropriate. |
| `w3lib-2.0.6` | Deprecated legacy | NCEP-origin | GRIB1-related support library. | Deprecated. | Identify dependencies and remove progressively. |

## Component priorities

The active development priority is:

1. `libmisc`;
2. `libax`;
3. `sigioBAM`, only as reviewed legacy code;
4. `w3lib-2.0.6`, only as temporary compatibility code.

## Core libraries

### libmisc

`libmisc` should be treated as the base utility library of the repository.

Expected characteristics:

- minimal dependencies;
- reusable routines;
- portable Fortran code;
- no direct dependence on model-specific formats;
- suitable as the first target for CMake and tests.

Initial inventory tasks:

- list all source files;
- identify public routines;
- identify routines used by `libax`;
- identify routines used by external applications;
- check whether each source file has license or authorship headers;
- separate generic routines from application-specific routines.

Expected classification rules:

| File type or routine type | Recommended classification |
| --- | --- |
| Generic string/date/file utilities | Keep in `libmisc`. |
| Scientific data access routines | Move or keep in `libax`. |
| BAM-specific routines | Review for `sigioBAM` or migration. |
| GRIB1/NCEP-specific routines | Review for removal or optional backend. |

### libax

`libax` should be treated as the scientific data access library.

Expected characteristics:

- support for meteorological and scientific data formats;
- possible dependencies on NetCDF or legacy GRIB-related code;
- clear dependency direction toward `libmisc`;
- gradual modernization toward a cleaner public API.

Initial inventory tasks:

- list all source files;
- identify supported formats;
- identify external dependencies;
- identify any remaining `w3lib-2.0.6` usage;
- document public routines and expected inputs/outputs;
- collect small sample files for future tests;
- separate format-specific code into clear backends when possible.

Expected format classification:

| Format or workflow | Recommended status |
| --- | --- |
| NetCDF | Active support. |
| GrADS | Active or legacy support, depending on current usage. |
| GRIB1 | Legacy support; isolate behind compatibility layer. |
| BAM-specific spectral files | Review whether it belongs in `libax` or `sigioBAM`. |

## Legacy libraries

### sigioBAM

`sigioBAM` is preserved because it may contain useful routines related to BAM spectral files.

Initial inventory tasks:

- list all source files;
- identify routines still used by active workflows;
- classify routines as keep, migrate, deprecate, or remove;
- identify routines that could move to `libax`;
- check authorship and license headers;
- document any dependency on `libmisc`, `libax`, or `w3lib-2.0.6`.

Recommended classification:

| Routine status | Meaning | Action |
| --- | --- | --- |
| Keep | Still useful and BAM-specific. | Preserve in `sigioBAM` while documenting. |
| Migrate | Useful beyond `sigioBAM`. | Move later to `libax` or `libmisc`. |
| Deprecate | Not recommended for new use. | Document replacement. |
| Remove | No longer used or duplicated. | Remove only after verification. |

### w3lib-2.0.6

`w3lib-2.0.6` is a deprecated NCEP-origin dependency.

Initial inventory tasks:

- identify whether it is still required;
- identify which source files call it;
- isolate any remaining dependency;
- preserve original license information while it remains in the repository;
- remove it only after confirming that no active workflow depends on it.

Recommended classification:

| Usage type | Action |
| --- | --- |
| Direct call from `libax` | Isolate behind GRIB1 compatibility layer. |
| Direct call from application code | Document and plan replacement. |
| Unused object/source file | Candidate for removal after verification. |
| License or notice file | Preserve while any `w3lib` code remains. |

## Intended dependency direction

The intended dependency direction is:

```text
applications
  |
  +-- libax
  |     |
  |     +-- libmisc
  |
  +-- libmisc
```

`libmisc` should not depend on `libax`, `sigioBAM`, or `w3lib-2.0.6`.

`libax` may depend on `libmisc` and external data access libraries.

Legacy dependencies should be isolated and documented.

## File-level inventory template

After the source tree is inspected, this document should include a file-level table like this:

| Path | Component | Language | Main routines | Dependencies | Status | License notes |
| --- | --- | --- | --- | --- | --- | --- |
| `path/to/file.F90` | `libmisc` | Fortran | To be identified | None expected | Active | Check header |
| `path/to/file.F90` | `libax` | Fortran | To be identified | NetCDF / legacy GRIB | Active | Check header |
| `path/to/file.F90` | `sigioBAM` | Fortran | To be identified | BAM-specific | Legacy | Check header |
| `path/to/file.f` | `w3lib-2.0.6` | Fortran | To be identified | NCEP legacy | Deprecated | Preserve original notice |

## Recommended local inspection commands

When working from a local clone, use:

```bash
git pull --ff-only
find . -maxdepth 3 -type d | sort
find . -maxdepth 4 -type f | sort
```

To identify source files:

```bash
find . -type f \( -name '*.f' -o -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c' -o -name '*.h' \) | sort
```

To identify possible `w3lib` usage:

```bash
grep -RIn "w3\|w3lib\|getgb\|putgb\|grib" . \
  --include='*.f' --include='*.F' --include='*.f90' --include='*.F90' \
  --include='*.c' --include='*.h'
```

To identify NetCDF usage:

```bash
grep -RIn "netcdf\|nf90_\|nf_" . \
  --include='*.f' --include='*.F' --include='*.f90' --include='*.F90' \
  --include='*.c' --include='*.h'
```

## Next step

The next development step should be a build-focused PR for `libmisc`.

Before that PR, this inventory should be updated with the actual file tree if a complete source listing becomes available.

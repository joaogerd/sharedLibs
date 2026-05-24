# Source inventory

This document records the initial source inventory for `sharedLibs`.

It is a planning document. The inventory should be refined after inspecting the actual source tree, build scripts, and applications that depend on these libraries.

## Component summary

| Component | Category | Current role | Development status | Recommended action |
| --- | --- | --- | --- | --- |
| `libmisc` | Core | General-purpose utility library. | Active target. | Consolidate, document, build, and test first. |
| `libax` | Core | Scientific and meteorological data access library. | Active target. | Document API, dependencies, and supported formats. |
| `sigioBAM` | Legacy | BAM spectral-file routines. | Maintained legacy. | Review routines and migrate useful code when appropriate. |
| `w3lib` | Deprecated legacy | NCEP-origin GRIB1-related library. | Deprecated. | Identify dependencies and remove progressively. |

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
- check whether each source file has license or authorship headers.

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
- identify any remaining `w3lib` usage;
- document public routines and expected inputs/outputs;
- collect small sample files for future tests.

## Legacy libraries

### sigioBAM

`sigioBAM` is preserved because it may contain useful routines related to BAM spectral files.

Initial inventory tasks:

- list all source files;
- identify routines still used by active workflows;
- classify routines as keep, migrate, deprecate, or remove;
- identify routines that could move to `libax`;
- check authorship and license headers.

### w3lib

`w3lib` is a deprecated NCEP-origin dependency.

Initial inventory tasks:

- identify whether it is still required;
- identify which source files call it;
- isolate any remaining dependency;
- preserve original license information while it remains in the repository;
- remove it only after confirming that no active workflow depends on it.

## Dependency direction

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

`libmisc` should not depend on `libax`, `sigioBAM`, or `w3lib`.

`libax` may depend on `libmisc` and external data access libraries.

Legacy dependencies should be isolated and documented.

## Next inventory step

The next step is to inspect the actual source tree and replace this planning inventory with a file-level table.

A future table should include:

```text
Path | Component | Language | Main routines | Dependencies | Status | License notes
```

This will support later build-system work and safer code modernization.

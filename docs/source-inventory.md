# Source inventory

This document records the current source inventory for `sharedLibs`.

The repository is now indexed for code search, and the inventory below reflects files confirmed from the current tree and from the existing Automake build metadata.

## Confirmed repository layout

The current source tree is organized mainly as:

```text
sharedLibs/
├── CPTECLIBS/
│   ├── libmisc/
│   ├── libax/
│   └── sigioBAM/
├── w3lib-2.0.6/
├── docs/
└── README.md
```

This confirms that the active CPTEC-related libraries are currently under `CPTECLIBS/`, while the NCEP-origin `w3lib-2.0.6` is located at the repository root.

## Confirmed repository components

| Component | Path | Category | Origin | Current role | Development status | Recommended action |
| --- | --- | --- | --- | --- | --- | --- |
| `libmisc` | `CPTECLIBS/libmisc/` | Core | Author-developed / CPTEC-related | General-purpose utility library. | Active target. | Consolidate, document, build, and test first. |
| `libax` | `CPTECLIBS/libax/` | Core | Author-developed / CPTEC-related | Scientific and meteorological data access library. | Active target. | Document API, dependencies, supported formats, and `w3lib` usage. |
| `sigioBAM` | `CPTECLIBS/sigioBAM/` | Legacy | BAM-related local code | Spectral-file routines for BAM workflows. | Maintained legacy. | Review routines and migrate useful code when appropriate. |
| `w3lib-2.0.6` | `w3lib-2.0.6/` | Deprecated legacy | NCEP-origin | GRIB1-related support library. | Deprecated. | Identify dependencies and remove progressively. |

## Build-system evidence

The three CPTEC-related libraries currently use an Autotools/Automake-style build structure.

Confirmed files include:

| Component | Confirmed build files |
| --- | --- |
| `libmisc` | `CPTECLIBS/libmisc/configure.ac`, `CPTECLIBS/libmisc/Makefile.am`, `CPTECLIBS/libmisc/src/Makefile.am`, generated `Makefile.in` files, `autogen.sh`, `aclocal.m4`, `m4/` macros. |
| `libax` | `CPTECLIBS/libax/configure.ac`, `CPTECLIBS/libax/Makefile.am`, `CPTECLIBS/libax/src/Makefile.am`, generated `Makefile.in` files, `autogen.sh`, `aclocal.m4`, `m4/` macros. |
| `sigioBAM` | `CPTECLIBS/sigioBAM/configure.ac`, `CPTECLIBS/sigioBAM/Makefile.am`, `CPTECLIBS/sigioBAM/src/Makefile.am`, generated `Makefile.in` files, `autogen.sh`, `aclocal.m4`, `m4/` macros. |

This means the first modernization step should not start by deleting the legacy build system. A safer approach is to document the current Autotools build, then add a modern CMake build in parallel, beginning with `libmisc`.

## Confirmed source files by component

### libmisc

The confirmed `libmisc` source list comes from `CPTECLIBS/libmisc/src/Makefile.am`.

The current library target is:

```make
lib_LTLIBRARIES = libmisc.la
```

The confirmed source files are:

| Source file | Component | Language | Initial classification | Notes |
| --- | --- | --- | --- | --- |
| `CPTECLIBS/libmisc/src/TypeKinds.f90` | `libmisc` | Fortran | Core utility | Type/kind definitions. |
| `CPTECLIBS/libmisc/src/m_inpak90.F90` | `libmisc` | Fortran | Core utility | Input/package-style utility; inspect API. |
| `CPTECLIBS/libmisc/src/EndianUtility.f90` | `libmisc` | Fortran | Core utility | Endianness helper. |
| `CPTECLIBS/libmisc/src/m_stdio.f90` | `libmisc` | Fortran | Core utility | Standard I/O helper. |
| `CPTECLIBS/libmisc/src/m_string.f90` | `libmisc` | Fortran | Core utility | String manipulation helper. |
| `CPTECLIBS/libmisc/src/m_time.f90` | `libmisc` | Fortran | Core utility | Time/date helper. |
| `CPTECLIBS/libmisc/src/m_msg.f90` | `libmisc` | Fortran | Core utility | Message/log helper. |
| `CPTECLIBS/libmisc/src/coord_compute.f90` | `libmisc` | Fortran | Review | Coordinate computation may overlap with `libax` or model/grid logic. |
| `CPTECLIBS/libmisc/src/BilinInterp.f90` | `libmisc` | Fortran | Review | Bilinear interpolation may remain as a generic numerical utility or move to a scientific utilities layer. |

Additional files found by code search and not yet classified:

| Source file | Initial action |
| --- | --- |
| `CPTECLIBS/libmisc/src/dateTimeMod.f90` | Inspect whether it is unused, legacy, or missing from the current build list. |
| `CPTECLIBS/libmisc/src/dateTimeMod_v2.f90` | Inspect relation with `m_time.f90`. |
| `CPTECLIBS/libmisc/src/sortingModule.f90` | Inspect whether it should be included in the active build. |
| `CPTECLIBS/libmisc/src/old/teste.f90` | Treat as old/test code; candidate for removal or migration to `tests/` after review. |

### libax

The confirmed `libax` source list comes from `CPTECLIBS/libax/src/Makefile.am`.

The current library target is:

```make
lib_LTLIBRARIES = libax.la
```

The confirmed source files are:

| Source file | Component | Language | Initial classification | Notes |
| --- | --- | --- | --- | --- |
| `CPTECLIBS/libax/src/coord_compute.f90` | `libax` | Fortran | Review | Duplicates or overlaps with `libmisc/src/coord_compute.f90`; inspect before modernization. |
| `CPTECLIBS/libax/src/accessGrib.F90` | `libax` | Fortran | Legacy backend | GRIB1/GRIB-related support; likely depends on `w3lib-2.0.6` or similar routines. |
| `CPTECLIBS/libax/src/accessNetcdf.f90` | `libax` | Fortran | Active backend | NetCDF access layer. |
| `CPTECLIBS/libax/src/m_GrADSfiles.F90` | `libax` | Fortran | Active or legacy backend | GrADS file support; inspect current use. |
| `CPTECLIBS/libax/src/fileAccess.f90` | `libax` | Fortran | Core interface | Likely the main file access abstraction. |

Initial design interpretation:

```text
libax
├── fileAccess.f90       # central interface candidate
├── accessNetcdf.f90     # NetCDF backend
├── m_GrADSfiles.F90     # GrADS backend
├── accessGrib.F90       # GRIB legacy backend
└── coord_compute.f90    # coordinate helper; possible duplication with libmisc
```

### sigioBAM

The confirmed `sigioBAM` source list comes from `CPTECLIBS/sigioBAM/src/Makefile.am`.

The current library target is:

```make
lib_LTLIBRARIES = libsigiobam.la
```

The confirmed source files are:

| Source file | Component | Language | Initial classification | Notes |
| --- | --- | --- | --- | --- |
| `CPTECLIBS/sigioBAM/src/coord_compute.F90` | `sigioBAM` | Fortran | Review | May duplicate coordinate logic from `libmisc` or `libax`. |
| `CPTECLIBS/sigioBAM/src/TypeKinds.f90` | `sigioBAM` | Fortran | Review | Duplicates type/kind logic from `libmisc`. |
| `CPTECLIBS/sigioBAM/src/EndianUtility.f90` | `sigioBAM` | Fortran | Review | Duplicates utility logic from `libmisc`. |
| `CPTECLIBS/sigioBAM/src/MiscMod.f90` | `sigioBAM` | Fortran | Legacy utility | Review for migration into `libmisc` if generic. |
| `CPTECLIBS/sigioBAM/src/ModConstants.f90` | `sigioBAM` | Fortran | BAM-specific or numerical constants | Inspect scope. |
| `CPTECLIBS/sigioBAM/src/LegendreTransform.f90` | `sigioBAM` | Fortran | Legacy numerical/BAM | Preserve until reviewed. |
| `CPTECLIBS/sigioBAM/src/Fourier.f90` | `sigioBAM` | Fortran | Legacy numerical/BAM | Preserve until reviewed. |
| `CPTECLIBS/sigioBAM/src/TransformTools.f90` | `sigioBAM` | Fortran | Legacy numerical/BAM | Preserve until reviewed. |
| `CPTECLIBS/sigioBAM/src/sigio_BAMMod.F90` | `sigioBAM` | Fortran | Main BAM spectral I/O | Preserve as legacy entry point. |

Additional files found by code search and not yet classified:

| Source file | Initial action |
| --- | --- |
| `CPTECLIBS/sigioBAM/src/Mod_LegendreTransform.f90` | Inspect relation with `LegendreTransform.f90`. |

### w3lib-2.0.6

`w3lib-2.0.6` contains many Fortran source files and GRIB1-related utilities.

Confirmed examples include:

| Source file | Initial classification |
| --- | --- |
| `w3lib-2.0.6/getgb1re.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/getgbeh.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/getgbmh.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/putgb.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/putgbex.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/putgben.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/putgbens.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/putgbn.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/skgb.f` | Deprecated GRIB1 support. |
| `w3lib-2.0.6/w3kind.f` | Deprecated NCEP support. |
| `w3lib-2.0.6/w3log.f` | Deprecated NCEP support. |
| `w3lib-2.0.6/w3movdat.f` | Deprecated date/time support. |
| `w3lib-2.0.6/w3locdat.f` | Deprecated date/time support. |
| `w3lib-2.0.6/w3utcdat.f` | Deprecated date/time support. |
| `w3lib-2.0.6/w3difdat.f` | Deprecated date/time support. |
| `w3lib-2.0.6/grib1.doc` | Documentation; preserve while GRIB1 support exists. |
| `w3lib-2.0.6/README` | Documentation; preserve while directory exists. |

`w3lib-2.0.6` should not be extended. The main task is to identify which routines are still called by `libax`, especially from `accessGrib.F90`, and isolate them behind a compatibility layer until they can be removed.

## Duplicated or overlapping code

The current tree contains several likely duplications:

| Pattern | Components involved | Action |
| --- | --- | --- |
| `coord_compute` appears in `libmisc`, `libax`, and `sigioBAM`. | `libmisc`, `libax`, `sigioBAM` | Compare implementations and decide whether one canonical implementation should live in `libmisc` or `libax`. |
| `TypeKinds.f90` appears in `libmisc` and `sigioBAM`. | `libmisc`, `sigioBAM` | Prefer one canonical definition in `libmisc` if compatible. |
| `EndianUtility.f90` appears in `libmisc` and `sigioBAM`. | `libmisc`, `sigioBAM` | Prefer one canonical implementation in `libmisc` if compatible. |
| Date/time utilities exist in `libmisc` and `w3lib-2.0.6`. | `libmisc`, `w3lib-2.0.6` | Prefer `libmisc`; avoid new `w3lib` date/time dependencies. |

## Repository hygiene observations

Code search revealed version-control artifacts inside the repository, especially `.svn/pristine/...` files under old CPTEC components.

These files should be treated as repository hygiene issues:

- they are not source files;
- they can confuse inventory and search results;
- they increase repository size;
- they should not be part of a clean Git repository unless there is a specific historical reason.

Recommended future action:

1. confirm that `.svn/` directories are not needed;
2. remove `.svn/` directories in a dedicated cleanup PR;
3. add `.svn/` to `.gitignore`;
4. avoid mixing cleanup with build-system modernization.

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

## Immediate technical implications

The inventory suggests the following order of work:

1. keep the current Autotools files for compatibility;
2. add a modern CMake build in parallel, starting only with `libmisc`;
3. validate `libmisc` sources independently;
4. only after that, add `libax` to the modern build;
5. treat `accessGrib.F90` as a legacy backend;
6. postpone `sigioBAM` modernization until duplications are understood;
7. remove `.svn/` directories in a separate cleanup PR.

## Recommended local inspection commands

When working from a local clone, use:

```bash
git pull --ff-only
find CPTECLIBS -maxdepth 4 -type f | sort
find w3lib-2.0.6 -maxdepth 2 -type f | sort
```

To identify source files:

```bash
find CPTECLIBS w3lib-2.0.6 -type f \
  \( -name '*.f' -o -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c' -o -name '*.h' \) | sort
```

To identify possible `w3lib` usage:

```bash
grep -RIn "w3\|w3lib\|getgb\|putgb\|grib" CPTECLIBS \
  --include='*.f' --include='*.F' --include='*.f90' --include='*.F90' \
  --include='*.c' --include='*.h'
```

To identify NetCDF usage:

```bash
grep -RIn "netcdf\|nf90_\|nf_" CPTECLIBS \
  --include='*.f' --include='*.F' --include='*.f90' --include='*.F90' \
  --include='*.c' --include='*.h'
```

## Next step

The next development step should be a build-focused PR for `libmisc`.

That PR should add a minimal CMake workflow for the confirmed `libmisc` source list while preserving the current Autotools files.

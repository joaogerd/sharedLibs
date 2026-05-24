# sharedLibs overview

`sharedLibs` is a collection of shared scientific libraries maintained to support meteorological, numerical modeling, data processing, and verification applications.

The repository is being reorganized around a smaller and clearer scope. The central components are `libmisc` and `libax`, while legacy components are preserved only when they still provide useful routines or compatibility with older workflows.

## Main purpose

The main purpose of this repository is to avoid duplicated utility code across scientific applications and to provide a reusable foundation for Fortran-based tools developed in CPTEC/INPE-related workflows.

The project should evolve toward:

- clear library boundaries;
- explicit documentation for each component;
- a reproducible build system;
- small examples showing how each library is used;
- tests for the most important routines;
- reduced dependence on external legacy libraries.

## Current library groups

### Core libraries

The core libraries are the parts of the repository that should receive active development.

- `libmisc`: general-purpose utility routines.
- `libax`: data access routines for scientific and meteorological file formats.

### Legacy libraries

Legacy libraries are preserved for compatibility or because they contain routines that may still be useful.

- `sigioBAM`: routines related to BAM spectral files.
- `w3lib`: NCEP-origin legacy library, planned for progressive removal.

## Development direction

The recommended development strategy is incremental.

First, the repository should document its current state and make the distinction between active and legacy components clear. Then, `libmisc` should become the first library with a clean build and test workflow. After that, `libax` should be reviewed and modernized, especially where it depends on legacy GRIB1-related code.

## Relationship with older CPTEC workflows

This repository may still be useful for applications such as BAM-related tools, SCANTEC, GSI utilities, pyBAM, and other local scientific programs. However, the goal is not only to preserve older code. The goal is to extract, organize, and improve reusable routines so that they can support newer workflows more safely.

## Relationship with NCEPLIBS

`sharedLibs` is not a mirror of NCEPLIBS and should not be documented as a replacement for NCEPLIBS.

External libraries from NCEP/NOAA should be treated as external dependencies or temporary compatibility components. The current long-term plan is to remove `w3lib` from the active dependency chain.

## Expected future layout

A possible target layout is:

```text
sharedLibs/
├── README.md
├── LICENSE
├── CMakeLists.txt
├── libs/
│   ├── libmisc/
│   ├── libax/
│   └── legacy/
│       └── sigioBAM/
├── docs/
├── examples/
└── tests/
```

The actual repository layout may still differ from this target while the reorganization is in progress.

# sharedLibs

`sharedLibs` is a collection of shared scientific libraries developed to support meteorological, numerical modeling, data processing, and verification applications.

The project is being reorganized around two main libraries:

- `libmisc`: general-purpose utilities used by scientific and operational codes.
- `libax`: data access utilities for meteorological and scientific file formats.

Most of the code in this repository is intended to be maintained as author-developed infrastructure for CPTEC/INPE-related applications and personal scientific software.

## Project direction

The long-term goal of `sharedLibs` is to provide a clean, documented, and maintainable set of reusable libraries for Fortran-based scientific applications.

The current development priority is:

1. consolidate and improve `libmisc`;
2. consolidate and modernize `libax`;
3. preserve useful legacy routines from `sigioBAM` while deciding what should be migrated or redesigned;
4. progressively remove the dependency on `w3lib`.

## Library overview

### libmisc

`libmisc` is the general utility library of the project.

It is intended to contain reusable routines that are independent of heavy external dependencies, such as:

- string manipulation;
- date and time handling;
- file and path checks;
- small numerical or logical utilities;
- common support routines used by other libraries and applications.

This library should remain simple, portable, and easy to test. It is the natural foundation for the other components of `sharedLibs`.

### libax

`libax` is the data access library of the project.

Its purpose is to provide routines for reading and handling scientific and meteorological datasets, including legacy and operational formats used in CPTEC/INPE workflows.

Historically, `libax` has been associated with access to formats such as:

- GrADS;
- GRIB1-based workflows;
- NetCDF.

The long-term goal is to evolve `libax` into a cleaner abstraction layer for scientific data access, hiding format-specific details behind a more consistent interface.

### sigioBAM

`sigioBAM` is currently kept as a legacy component.

It contains routines related to spectral files from the BAM model and still preserves useful code that may be improved, reused, or migrated in the future.

Although `sigioBAM` may eventually stop existing as an independent library, it should not be removed until its useful routines have been reviewed and either migrated or explicitly deprecated.

### w3lib

`w3lib` is the only NCEP-origin library currently associated with this repository.

It is considered deprecated in the context of `sharedLibs` and should be progressively removed from the project. New development should not depend on it unless strictly necessary for temporary compatibility with legacy workflows.

## Suggested repository organization

The repository is expected to move toward a structure similar to:

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
│   ├── overview.md
│   ├── libmisc.md
│   ├── libax.md
│   ├── sigioBAM.md
│   └── migration-from-w3lib.md
├── examples/
│   ├── libmisc/
│   └── libax/
└── tests/
    ├── libmisc/
    └── libax/
```

This structure is a target organization. The current code may still reflect an older layout.

## Build system

A unified build system is planned for the repository.

The preferred direction is to provide an out-of-source CMake workflow, for example:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/sharedLibs
cmake --build build -j
cmake --install build
```

At this stage, users should inspect the individual library directories for available build scripts or legacy compilation instructions.

## Development priorities

The recommended roadmap is:

1. document the current purpose of each library;
2. define the license and authorship status of the repository;
3. create a minimal build workflow for `libmisc`;
4. add small tests for `libmisc` routines;
5. document the current API of `libax`;
6. identify where `libax` still depends on `w3lib`;
7. isolate or replace `w3lib` usage;
8. review `sigioBAM` and migrate useful routines when appropriate.

## Relationship with NCEPLIBS

Earlier versions of this repository included broader documentation about NCEPLIBS.

The current direction is different: `sharedLibs` is not intended to be a mirror or replacement for NCEPLIBS. External libraries from NCEP/NOAA should be treated as external dependencies or legacy compatibility components when needed.

## License

The license status of the repository should be reviewed and made explicit in a top-level `LICENSE` file.

Because this repository may contain code with different origins, each component should be checked before assigning or changing its license.

## Status

This repository is under reorganization.

The main active targets are `libmisc` and `libax`. Legacy components are preserved only while they remain useful for compatibility, migration, or future refactoring.

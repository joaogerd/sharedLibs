# Build guide

This document describes the current and planned build workflows for `sharedLibs`.

## Current status

The legacy CPTEC libraries still contain Autotools/Automake files.

The modernization path is to add a CMake workflow in parallel, starting with `libmisc`, while preserving the existing Autotools files until the new workflow is validated.

## CMake build for libmisc

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

To install into a local prefix:

```bash
cmake --install build --prefix "$HOME/sharedLibs"
```

## Scope

The initial CMake workflow builds only `libmisc`.

`libax`, `sigioBAM`, and `w3lib-2.0.6` should not be added to the modern build until their dependencies and legacy interactions are fully reviewed.

## Build policy

- Keep the legacy Autotools files while CMake is introduced.
- Start with the smallest stable target: `libmisc`.
- Add `libax` only after `libmisc` builds reliably.
- Treat GRIB1 and `w3lib-2.0.6` support as legacy compatibility.
- Do not modernize `sigioBAM` before duplicated utility code is reviewed.

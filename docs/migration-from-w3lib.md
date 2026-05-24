# Migration from w3lib

`w3lib` is the only NCEP-origin library currently associated with `sharedLibs`.

In the current direction of the project, `w3lib` is treated as a legacy dependency. New development should avoid depending on it unless temporary compatibility with older workflows is required.

## Current status

`w3lib` is deprecated in the context of `sharedLibs`.

This does not mean it must be removed immediately. It means the project should not expand its use or treat it as a core component.

## Why remove it gradually

The reasons for gradually removing `w3lib` are:

- it is not an author-developed component of this repository;
- it increases the legacy footprint of the project;
- it makes the identity of `sharedLibs` less clear;
- it can complicate future build and licensing decisions.

## Migration strategy

The migration should be incremental.

Recommended steps:

1. identify all files that call routines from `w3lib`;
2. identify which libraries depend on these calls, especially `libax`;
3. classify each use as required, replaceable, or obsolete;
4. isolate the remaining calls behind a small compatibility layer;
5. replace direct calls where possible;
6. remove unused code only after confirming that no active workflow depends on it.

## Relationship with libax

If `libax` still depends on `w3lib`, that dependency should be documented explicitly.

The preferred intermediate solution is to isolate GRIB1-specific routines so that the rest of `libax` does not depend directly on `w3lib`.

A possible organization is:

```text
libax core routines
  NetCDF support
  GrADS support
  legacy GRIB1 support through a compatibility layer
```

This keeps the legacy dependency contained while allowing the rest of `libax` to evolve.

## Replacement options

Replacement options should be evaluated only after inspecting the current code.

Possible paths include:

- replacing legacy GRIB1 support with a maintained external library;
- preserving only the minimal routines required for older workflows;
- dropping GRIB1 support if no active workflow depends on it;
- moving GRIB1 support to an optional backend.

## Documentation needed

The next documentation step is to create an inventory table:

```text
Source file | w3lib routine | Called by | Required? | Migration action
```

This table should guide the removal process and reduce the risk of breaking working applications.

## Policy for new code

New code should not call `w3lib` directly.

If temporary use is unavoidable, it should be documented as legacy compatibility and isolated from the main API.

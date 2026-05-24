# libmisc

`libmisc` is the general-purpose utility library of `sharedLibs`.

It should provide small, reusable routines that can be shared by scientific and operational applications without introducing heavy external dependencies.

## Role in the project

`libmisc` is expected to be the most basic and portable library in the repository.

It should be usable by other libraries, including `libax`, and by external applications that need common support routines.

Because of this role, `libmisc` should remain simple, well documented, and easy to test.

## Expected scope

The library may include routines for:

- string manipulation;
- date and time handling;
- file and path checks;
- filename formatting;
- small numerical utilities;
- status and error handling helpers;
- common constants or simple configuration helpers.

## Design principles

`libmisc` should follow a few basic principles:

1. keep dependencies minimal;
2. avoid application-specific assumptions;
3. prefer small routines with clear responsibilities;
4. document each public routine;
5. add tests for routines that are reused by other components;
6. preserve backward compatibility when possible.

## Dependency policy

`libmisc` should not depend on NetCDF, GRIB libraries, MPI, or model-specific code.

If a routine requires a heavy external library, it probably belongs in another component, not in `libmisc`.

## Build priority

`libmisc` should be the first library to receive a clean build workflow.

A future CMake target could be named:

```cmake
sharedlibs_misc
```

or simply:

```cmake
misc
```

The final target name should be chosen after inspecting the existing code and avoiding conflicts with legacy build scripts.

## Testing priority

The first tests should focus on deterministic routines such as:

- string trimming and formatting;
- date conversions;
- file existence checks;
- path manipulation;
- simple numerical helpers.

These tests can be small, but they are important because `libmisc` may become a dependency of other libraries.

## Modernization notes

When modernizing `libmisc`, avoid rewriting everything at once.

A safer approach is:

1. identify the public routines currently used by other applications;
2. document their expected behavior;
3. add small tests around the existing behavior;
4. refactor internally only after behavior is protected by tests.

This allows the library to improve without breaking older codes unexpectedly.

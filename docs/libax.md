# libax

`libax` is the data access library of `sharedLibs`.

Its purpose is to provide reusable routines for reading, inspecting, and handling scientific and meteorological datasets used in CPTEC/INPE-related workflows.

## Role in the project

`libax` should become the main abstraction layer for data access.

Instead of each application implementing its own file-reading logic, `libax` should concentrate common routines for working with meteorological and scientific data formats.

Historically, this includes workflows involving:

- GrADS files;
- GRIB1-based files;
- NetCDF files;
- legacy files used by CPTEC applications.

## Long-term goal

The long-term goal is to provide a cleaner and more consistent interface for data access.

A future interface should help applications answer questions such as:

- what is the file format?
- which variables are available?
- what are the dimensions of the dataset?
- what is the grid structure?
- how can a 2D or 3D field be read?
- how can format-specific details be hidden from the calling application?

## Expected responsibilities

`libax` may include routines for:

- opening scientific data files;
- identifying file formats;
- reading metadata;
- reading fields;
- handling grids and dimensions;
- converting or exposing data in a common internal representation;
- supporting legacy CPTEC workflows while they are still needed.

## Dependency policy

Unlike `libmisc`, `libax` may depend on external scientific I/O libraries.

Possible dependencies include:

- NetCDF C/Fortran libraries;
- GRIB-related libraries, if still required;
- other format-specific dependencies.

These dependencies should be documented clearly and isolated as much as possible.

## Relationship with libmisc

`libax` may use `libmisc` for general utilities such as string handling, date processing, file checks, and common error handling.

However, `libmisc` should not depend on `libax`.

The intended dependency direction is:

```text
libax -> libmisc
```

not the opposite.

## Relationship with w3lib

Any dependency on `w3lib` should be considered temporary.

If `libax` currently uses `w3lib` for GRIB1-related workflows, this usage should be documented and isolated. The long-term goal is to either replace it or keep it behind a narrow compatibility layer until it can be removed.

## Modernization strategy

`libax` should be modernized carefully because it may touch older operational workflows.

A safe strategy is:

1. document the current public routines;
2. identify supported formats and required dependencies;
3. identify which applications use each routine;
4. isolate legacy GRIB1 and `w3lib` usage;
5. add small example programs for each supported format;
6. add tests using small sample files;
7. gradually introduce a cleaner public API.

## Possible future API direction

A future high-level API could provide calls similar to:

```fortran
call ax_open(filename, handler)
call ax_get_variables(handler, variables)
call ax_get_dimensions(handler, dimensions)
call ax_read_field(handler, variable_name, field)
call ax_close(handler)
```

This is only a conceptual direction. The final API should be designed after reviewing the current implementation and the applications that already depend on it.

## Documentation needs

The next documentation steps for `libax` are:

- list all source files;
- list public routines;
- identify supported formats;
- document dependencies;
- document examples of use;
- identify routines that should be deprecated or migrated.

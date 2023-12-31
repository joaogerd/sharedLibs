# sharedLibs


This is a collection of commonly shared libraries known as sharedLibs, which are required for several CPTEC applications. There are two main types of libraries: those developed by CPTEC (referred to as CPTECLIBS) and those developed by NCEP (referred to as NCEPLIBS).

The NCEPLIBS were forked from the [NOAA-EMC](https://github.com/NOAA-EMC) GitHub repository, and all documentation can be accessed there.

Below is a brief description of CPTECLIBS and NCEPLIBS (NCEPLIBS description is from the [NOAA-EMC](https://github.com/NOAA-EMC) GitHub readme).


## CPTECLIBS

### Introduction

This is a collection of libraries commonly known as CPTECLIBS that are
required for several CPTEC applications e.g. GSI, SCANTEC, pyBAM etc. For
general information about CPTECLIBS, the user is referred to the
[Wiki](https://github.com/joaogerd/CPTECLIBS/wiki).


## NCEPLIBS

### Introduction

This is a collection of libraries commonly known as NCEPLIBS that are
required for several NCEP applications e.g. UFS, GSI, UPP, etc. For
general information about NCEPLIBS, the user is referred to the
[Wiki](https://github.com/NOAA-EMC/NCEPLIBS/wiki).

For full documentation of each package see
https://noaa-emc.github.io/NCEPLIBS/.

For installing NCEPLIBS on NOAA HPC systems, see the
[hpc-stack](https://github.com/NOAA-EMC/hpc-stack) project. To install
the third-party libraries needed by NCEPLIBS, see
[NCEPLIBS-external](https://github.com/NOAA-EMC/NCEPLIBS-external).

NCEPLIBS is comprised of the following libraries:

 NCEP library         | Description
----------------------|---------------------------------------------------
[bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)             | binary I/O for the NCEP models
[bufr](https://github.com/NOAA-EMC/NCEPLIBS-bufr)               | working with the BUFR format
[g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)                   | coding/decoding GRIB2 messages
[g2c](https://github.com/NOAA-EMC/NCEPLIBS-g2c)                 | C decoder/encoder routines for GRIB edition 2
[g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl)           | Utilities for GRIB2 templates
[gfsio](https://github.com/NOAA-EMC/NCEPLIBS-gfsio)             | convert GFS Gaussian output into grib output
[ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)                   | general interpolation library 
[landsfcutil](https://github.com/NOAA-EMC/NCEPLIBS-landsfcutil) | initializing land-surface states
[ncio](https://github.com/NOAA-EMC/NCEPLIBS-ncio)               | read the GFS forecast files
[nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio)           | I/O for the NCEP models using NEMS
[nemsiogfs](https://github.com/NOAA-EMC/NCEPLIBS-nemsiogfs)     | I/O for the NEMS-GFS model
[sfcio](https://github.com/NOAA-EMC/NCEPLIBS-sfcio)             | surface files I/O
[sigio](https://github.com/NOAA-EMC/NCEPLIBS-sigio)             | I/O on the sigma restart file of the global spectral model
[sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)                   | spectral transform library 
[w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)             | decoder/encoder routines for GRIB edition 1
[w3nco](https://github.com/NOAA-EMC/NCEPLIBS-w3nco)             | decoder/encoder routines for GRIB edition 1
[wrf_io](https://github.com/NOAA-EMC/NCEPLIBS-wrf_io)           | lightweight WRF-IO API library for Unified Post Processor (UPP)
[ip2](https://github.com/NOAA-EMC/NCEPLIBS-ip2)                 | general interpolation library 2

The following projects are also related to the NCEPLIBS project:
 - [NCEPLIBS-grib_util](https://github.com/NOAA-EMC/NCEPLIBS-grib_util)
 - [NCEPLIBS-prod_util](https://github.com/NOAA-EMC/NCEPLIBS-prod_util)
 - [UFS_UTILS](https://github.com/NOAA-EMC/UFS_UTILS)

### Required Software

A supported C/C++ and Fortran compiler (see table below).  Other
versions may work, in particular if close to the versions listed
below.

| Compiler vendor | Supported (tested) versions                                |
|-----------------|------------------------------------------------------------|
| Intel           | 18.0.3.222, 18.0.5.274, 19.0.2.187, 19.0.5.281, 19.1.0.166 |
| GNU             | 8.3.0, 9.X.0                                               |

A supported MPI library (see table below).  Other versions may work,
in particular if close to the versions listed below.

| MPI library     | Supported (tested) versions                                |
|-----------------|------------------------------------------------------------|
| MPICH           | 3.3.1, 3.3.2                                               |
| MVAPICH2        | 2.3.3                                                      |
| Open MPI        | 3.1.5, 4.0.2                                               |
| Intel MPI       | 2018.0.4, 2019.6.154, 2020.0.166                           |
| SGI MPT         | 2.19                                                       |

A software stack consisting of the following third-party libraries
(TPL), compiled with the same compiler and MPI library (where
applicable). Other versions may work, in particular if close to the
versions listed below. To install the third-party libraries needed by
NCEPLIBS, see
[NCEPLIBS-external](https://github.com/NOAA-EMC/NCEPLIBS-external).

| Library         | Supported (tested) versions                                |
|-----------------|------------------------------------------------------------|
| CMake           | cmake-3.16.3                                               |
| zlib            | zlib-1.2.11                                                |
| libpng          | libpng-1.6.35                                              |
| libjpeg         | jpeg-9.1                                                   |
| Jasper          | jasper-2.0.16                                              |
| HDF5            | hdf5-1.10.4                                                |
| NetCDF          | netcdf-c-4.7.3, netcdf-fortran-4.5.2                       |

The [NCEPLIBS-external](https://github.com/NOAA-EMC/NCEPLIBS-external)
project provides an convenient way to install the required TPL on
supported platforms such as generic Linux/macOS systems.

### Obtaining, Building, Using NCEPLIBS

#### Obtaining the code

Clone the repository:

`git clone https://github.com/NOAA-EMC/NCEPLIBS`

#### Building the libraries

`CMake` employs an out-of-source build.  Create a directory for
configuring the build and cd into it:

`mkdir -p build && cd build`

Set the compilers, if needed, to match those being used for compiling
the TPL listed above: `FC`, `CC`, and `CXX` environment variables can
be used to point to the desired Fortran, C, and C++ compilers.

Execute `cmake` from inside your build directory.

`cmake -DCMAKE_INSTALL_PREFIX=<nceplibs-prefix> ../NCEPLIBS`

If the TPL are not located in a path recognized by `cmake`
e.g. `/usr/local`, it may be necessary to provide the appropriate
environment variables e.g. `<package_ROOT>` or `CMAKE_PREFIX_PATH` so
that `cmake` is able to locate these dependencies.

The installation prefix for NCEPLIBS is provided by the `cmake`
command-line argument `-DCMAKE_INSTALL_PREFIX=<nceplibs-prefix>`

To build and install:

```
make -j<x>
```

This will build and install all the NCEPLIBS and place them into the
installation prefix `<nceplibs-prefix>`.

Additional `cmake` command-line options while configuring NCEPLIBS:

`-DFLAT=ON|OFF` - determines the installation tree heirarchy

`-DDEPLOY=ON|OFF` - prepares and deploys Lua modulefiles for each of
the NCEPLIBS libraries.  To deploy modulefiles:

```
make deploy
```

[LMod](https://lmod.readthedocs.io/en/latest/) Environment module
system is required to load/unload these modules.

Some HPC systems use Tcl based Environment modules.  NCEPLIBS provides
a means to deploy Tcl module files and can be exercised as follows:

`-DTCLMOD=ON` - enables option to change Lua modulefiles to Tcl modulefiles for each of the NCEPLIBS libraries.  To make the change:
```
make tcl
```

On machines that do not provide access to Github, a `DOWNLOAD_ONLY`
and `USE_LOCAL` option is provided.  To download the NCEPLIBS without
building them use the cmake command line option `-DDOWNLOAD_ONLY=ON`.
This will download all the NCEPLIBS in a local `download` folder.
`make tarball` will tar and zip this folder into `nceplibs.tgz` which
can then be transferred to the machine without GitHub access.  Untar
the tarball, and use the cmake command line option `USE_LOCAL` as
`-DUSE_LOCAL=ON`.  This will use the source code from the download
directory.

#### Usage

`NCEPLIBS` can be used in any application that uses `cmake` to
configure and build by adding `-DCMAKE_PREFIX_PATH=<nceplibs-prefix>`
to the cmake command line during configuration.

### Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

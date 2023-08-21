# The lib_LIBRARIES macro tells the autotools the name and type of library we
# want to build - in this case _not_ a shared library. A shared library
# would end with '.la' instead of '.a'.
lib_LTLIBRARIES = libsigiobam.la
# This specially named variable contains a space delimited list of the source code
# files that must be compiled in order to create the library. The variable name
# is in the format:
# <library_name>_<library_type>_SOURCES
# Where:
# <library_name> = the name of the library specified in lib_LIBRARIES
# <library_type> = either 'a' for non-shared library or 'la' for shared.
libsigiobam_la_SOURCES = coord_compute.F90 TypeKinds.f90 EndianUtility.f90 MiscMod.f90 ModConstants.f90 LegendreTransform.f90 Fourier.f90 TransformTools.f90 sigio_BAMMod.F90
# Header files that must be installed are specified by the HEADERS family
# of variables.
include_HEADERS = coord_compute.mod typekinds.mod endianutility.mod miscmod.mod modconstants.mod legendretransform.mod fourier.mod transformtools.mod sigio_bammod.mod

libsigiobam_la_LDFLAGS=-version-info 1:0:0

CLEANFILES = *.$(FC_MODEXT)

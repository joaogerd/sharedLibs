#Automake uses libtool to build libraries declared with the LTLIBRARIES primary. 
#Each _LTLIBRARIES variable is a list of libtool libraries to build. 
lib_LTLIBRARIES = libax.la
# This specially named variable contains a space delimited list of the source code
# files that must be compiled in order to create the library. The variable name
# is in the format:
# <library_name>_<library_type>_SOURCES
# Where:
# <library_name> = the name of the library specified in lib_LIBRARIES
# <library_type> = either 'a' for non-shared library or 'la' for shared.
libax_la_SOURCES = coord_compute.f90 accessGrib.F90 accessNetcdf.f90 m_GrADSfiles.F90 fileAccess.f90
# Header files that must be installed are specified by the HEADERS family
# of variables.
#
include Makefile-headers
include Makefile-deps
#include_HEADERS = coord_compute.$(FC_MODEXT) accessgrib.$(FC_MODEXT) accessnetcdf.$(FC_MODEXT) m_gradsfiles.$(FC_MODEXT) fileaccess.$(FC_MODEXT)

libax_la_LDFLAGS=-version-info 1:0:0

CLEANFILES = *.$(FC_MODEXT)

#include Makefile.deps

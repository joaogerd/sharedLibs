cmake_minimum_required(VERSION 3.18)

project(sharedLibs VERSION 0.1.0 LANGUAGES Fortran)

option(SHAREDLIBS_BUILD_LIBMISC "Build the libmisc library" ON)

if(SHAREDLIBS_BUILD_LIBMISC)
  add_subdirectory(CPTECLIBS/libmisc)
endif()

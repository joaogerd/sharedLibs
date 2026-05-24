add_library(sharedlibs_misc)

target_sources(sharedlibs_misc
  PRIVATE
    src/TypeKinds.f90
    src/m_inpak90.F90
    src/EndianUtility.f90
    src/m_stdio.f90
    src/m_string.f90
    src/m_time.f90
    src/m_msg.f90
    src/coord_compute.f90
    src/BilinInterp.f90
)

set_target_properties(sharedlibs_misc PROPERTIES
  OUTPUT_NAME misc
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

target_include_directories(sharedlibs_misc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:include/sharedLibs/libmisc>
)

install(TARGETS sharedlibs_misc
  EXPORT sharedLibsTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
  DESTINATION include/sharedLibs/libmisc
  FILES_MATCHING PATTERN "*.mod"
)

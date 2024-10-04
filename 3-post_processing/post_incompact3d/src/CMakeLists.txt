include_directories(${2decomp_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src)	               

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(post_incompact3d
                    ../../../4-common/src/derive.f90
                    ../../../4-common/src/derive_y1d.f90
                    ../../../4-common/src/extra_tools.f90
                    ../../../4-common/src/filters.f90
                    ../../../4-common/src/ibm.f90
                    ../../../4-common/src/implicit.f90
                    ../../../4-common/src/module_param.f90
                    ../../../4-common/src/parameters.f90
                    post_incompact3d.f90
                    post_processing_mod.f90
                    ../../../4-common/src/schemes.f90
                    statistics_subs.f90
                    ../../../4-common/src/tools.f90
                    ../../../4-common/src/variables.f90
                    ../../../4-common/src/visu.f90)

target_link_libraries(post_incompact3d PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(post_incompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(ADIOS2_FOUND)
  target_link_libraries(post_incompact3d PRIVATE adios2::fortran_mpi)
endif()

install(TARGETS post_incompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

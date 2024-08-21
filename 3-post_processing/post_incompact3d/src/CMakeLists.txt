include_directories(${CMAKE_SOURCE_DIR}/decomp2d
	            ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(post_incompact3d
                    ../../4-incompact3d_common/derive.f90
                    ../../4-incompact3d_common/derive_y1d.f90
                    ../../4-incompact3d_common/extra_tools.f90
                    ../../4-incompact3d_common/filters.f90
                    ../../4-incompact3d_common/ibm.f90
                    ../../4-incompact3d_common/implicit.f90
                    ../../4-incompact3d_common/module_param.f90
                    ../../4-incompact3d_common/parameters.f90
                    post.f90
                    post_processing_mod.f90
                    ../../4-incompact3d_common/schemes.f90
                    statistics_subs.f90
                    tools_pp.f90
                    ../../4-incompact3d_common/variables.f90
                    ../../4-incompact3d_common/visu.f90)
                    
target_include_directories(post_incompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_link_libraries(post_incompact3d PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(post_incompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
if(ADIOS2_FOUND)
  #target_link_libraries(post_incompact3d ${ADIOS_LIBRARIES})
  target_link_libraries(post_incompact3d PRIVATE adios2::fortran_mpi)
endif()

install(TARGETS post_incompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

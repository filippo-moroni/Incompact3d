include_directories(${CMAKE_SOURCE_DIR}/decomp2d
	            ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(post_incompact3d
                    BC-dbg-schemes.f90
                    derive.f90
                    filters.f90
                    ibm.f90
                    module_param.f90
                    parameters.f90
                    post.f90
                    post_processing_mod.f90
                    schemes.f90               
                    statistics_sub.f90
                    tools.f90               
                    variables.f90
                    visu.f90)
                    
#add_executable(post_incompact3d ${files_post_incompact3d})
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

include_directories(../../4-common/decomp2d
	                  ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d 
                    BC-Channel-flow.f90                
                    BC-Temporal-TBL.f90
                    case.f90
                    ../../4-common/src/derive.f90
                    ../../4-common/src/derive_y1d.f90
                    ../../4-common/src/extra_tools.f90
                    ../../4-common/src/filters.f90
                    genepsi3d.f90
                    ../../4-common/src/ibm.f90
                    ../../4-common/src/implicit.f90
                    les_models.f90
                    ../../4-common/src/module_param.f90
                    navier.f90
                    ../../4-common/src/parameters.f90
                    poisson.f90
                    ../../4-common/src/schemes.f90
                    time_integrators.f90
                    tools.f90
                    transeq.f90                 
                    ../../4-common/src/variables.f90
                    ../../4-common/src/visu.f90
                    xcompact3d.f90)

target_include_directories(xcompact3d PRIVATE ${PROJECT_BINARY_DIR}/decomp2d)
target_link_libraries(xcompact3d PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(ADIOS2_FOUND)
  target_link_libraries(xcompact3d PRIVATE adios2::fortran_mpi)
endif()

install(TARGETS xcompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

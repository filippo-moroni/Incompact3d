include_directories(../../4-common/decomp2d
	            ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d 
                    BC-Channel-flow.f90                
                    BC-Temporal-TBL.f90
                    case.f90
                    constants.f90
                    ../../4-common/derive.f90
                    ../../4-common/derive_y1d.f90
                    ../../4-common/extra_tools.f90
                    ../../4-common/filters.f90
                    genepsi3d.f90
                    ../../4-common/ibm.f90
                    ../../4-common/implicit.f90
                    les_models.f90
                    ../../4-common/module_param.f90
                    navier.f90
                    ../../4-common/parameters.f90
                    poisson.f90
                    ../../4-common/schemes.f90
                    time_integrators.f90
                    tools.f90
                    transeq.f90                 
                    ../../4-common/variables.f90
                    ../../4-common/visu.f90
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

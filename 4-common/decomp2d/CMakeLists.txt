file(GLOB files_decomp decomp_2d.f90
                       glassman.f90
		                   io.f90)

# Generic FFT
if(${FFT_Choice} MATCHES "generic")
  file(GLOB files_fft fft_generic.f90)
endif(${FFT_Choice} MATCHES "generic")

# FFTW3
if(${FFT_Choice} MATCHES "fftw3")
  file(GLOB files_fft fft_fftw3.f90)
endif(${FFT_Choice} MATCHES "fftw3")

# FFT MKL (Intel)
if(${FFT_Choice} MATCHES "mkl")
  file(GLOB files_fft fft_mkl.f90)
endif(${FFT_Choice} MATCHES "mkl")

set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})

# Add the directory where the FFTW3 file 'fftw3.f' is located
if(${FFT_Choice} MATCHES "fftw3")
  target_include_directories(decomp2d PRIVATE ../../5-external/fftw-3.3.10/include)

  #target_include_directories(decomp2d PRIVATE ../../5-external/fftw-3.3.10)

  # Add the FFTW3 library directory
  target_link_directories(decomp2d PRIVATE ../../5-external/fftw-3.3.10/lib)

  # Link against the FFTW3 library
  target_link_libraries(decomp2d PRIVATE fftw3)
  
  endif(${FFT_Choice} MATCHES "fftw3")

if (MPI_FOUND)
  target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND) 

if(ADIOS2_FOUND)
  target_link_libraries(decomp2d PRIVATE adios2::fortran_mpi)
endif()

install(TARGETS decomp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
)

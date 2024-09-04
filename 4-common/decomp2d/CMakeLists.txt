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

# Create the static library for decomp2d
add_library(decomp2d STATIC ${SRCFILES})

if(${FFT_Choice} MATCHES "fftw3")
  
  # Add the directory where the FFTW3 file 'fftw3.f' is located (header)
  target_include_directories(decomp2d PRIVATE ../../5-external/fftw-3.3.10/include)

  # Add the FFTW3 directory where 'libfftw3.so' is located
  target_link_libraries(decomp2d PRIVATE ../../5-external/fftw-3.3.10/lib/libfftw3.so)
  #target_link_libraries(decomp2d PRIVATE ../../5-external/fftw-3.3.10/lib/libfftw3_threads.so.3.6.10)

endif(${FFT_Choice} MATCHES "fftw3")

#!---------------------------------------------------------------!

# Check if MKLROOT is defined
if(NOT DEFINED ENV{MKLROOT})
  message(FATAL_ERROR "MKLROOT environment variable is not set.")
endif()

# Add MKL FFT dependencies manually
if(${FFT_Choice} MATCHES "mkl")
  set(MKL_ROOT $ENV{MKLROOT})

  # Include the MKL headers
  target_include_directories(decomp2d PRIVATE ${MKL_ROOT}/include)

  # Link MKL libraries
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_link_libraries(decomp2d PRIVATE
      ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.a
      ${MKL_ROOT}/lib/intel64/libmkl_sequential.a
      ${MKL_ROOT}/lib/intel64/libmkl_core.a
      -lpthread -lm -ldl
    )
  else()
    message(FATAL_ERROR "MKL linking only supported with Intel compilers in this example.")
  endif()
endif()

#!---------------------------------------------------------------!

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

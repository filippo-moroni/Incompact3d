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

# MKL static linking
if(${FFT_Choice} MATCHES "mkl")

  # Add MKL include directory
  target_include_directories(decomp2d PRIVATE ${MKLROOT}/include)

  # Link the necessary static MKL libraries
  target_link_libraries(decomp2d PRIVATE ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a)
  target_link_libraries(decomp2d PRIVATE ${MKLROOT}/lib/intel64/libmkl_sequential.a)
  target_link_libraries(decomp2d PRIVATE ${MKLROOT}/lib/intel64/libmkl_core.a)
  target_link_libraries(decomp2d PRIVATE -lpthread -lm -ldl)
  
endif(${FFT_Choice} MATCHES "mkl")


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

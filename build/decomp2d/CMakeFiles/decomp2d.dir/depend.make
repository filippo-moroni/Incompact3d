# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o: \
 ../decomp2d/alloc.inc \
 ../decomp2d/factor.inc \
 ../decomp2d/halo.inc \
 ../decomp2d/halo_common.inc \
 ../decomp2d/transpose_x_to_y.inc \
 ../decomp2d/transpose_y_to_x.inc \
 ../decomp2d/transpose_y_to_z.inc \
 ../decomp2d/transpose_z_to_y.inc
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o: /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi/mpi.mod
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o.provides.build: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod decomp2d/decomp_2d.mod decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp GNU
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/build: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o: \
 ../decomp2d/fft_common.inc \
 ../decomp2d/fft_common_3d.inc
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o: decomp2d/CMakeFiles/decomp2d.dir/glassman.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o.provides.build: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_fft.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_fft.mod.stamp: decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod decomp2d/decomp_2d_fft.mod decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_fft.mod.stamp GNU
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/build: decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o.provides.build: decomp2d/CMakeFiles/decomp2d.dir/glassman.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/glassman.mod.stamp: decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod decomp2d/glassman.mod decomp2d/CMakeFiles/decomp2d.dir/glassman.mod.stamp GNU
decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/build: decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o: \
 ../decomp2d/io_read_inflow.f90 \
 ../decomp2d/io_read_one.inc \
 ../decomp2d/io_read_var.inc \
 ../decomp2d/io_write_every.inc \
 ../decomp2d/io_write_one.inc \
 ../decomp2d/io_write_outflow.f90 \
 ../decomp2d/io_write_plane.inc \
 ../decomp2d/io_write_var.inc
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o: /usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi/mpi.mod
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o.provides.build: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_io.mod.stamp
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_io.mod.stamp: decomp2d/CMakeFiles/decomp2d.dir/io.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod decomp2d/decomp_2d_io.mod decomp2d/CMakeFiles/decomp2d.dir/decomp_2d_io.mod.stamp GNU
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch decomp2d/CMakeFiles/decomp2d.dir/io.f90.o.provides.build
decomp2d/CMakeFiles/decomp2d.dir/build: decomp2d/CMakeFiles/decomp2d.dir/io.f90.o.provides.build

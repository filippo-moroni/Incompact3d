# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/n286654/Desktop/Incompact3d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/n286654/Desktop/Incompact3d/build

# Include any dependencies generated for this target.
include decomp2d/CMakeFiles/decomp2d.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include decomp2d/CMakeFiles/decomp2d.dir/compiler_depend.make

# Include the progress variables for this target.
include decomp2d/CMakeFiles/decomp2d.dir/progress.make

# Include the compile flags for this target's objects.
include decomp2d/CMakeFiles/decomp2d.dir/flags.make

decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o: decomp2d/CMakeFiles/decomp2d.dir/flags.make
decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o: ../decomp2d/decomp_2d.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/n286654/Desktop/Incompact3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/n286654/Desktop/Incompact3d/decomp2d/decomp_2d.f90 -o CMakeFiles/decomp2d.dir/decomp_2d.f90.o

decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/decomp2d.dir/decomp_2d.f90.i"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/n286654/Desktop/Incompact3d/decomp2d/decomp_2d.f90 > CMakeFiles/decomp2d.dir/decomp_2d.f90.i

decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/decomp2d.dir/decomp_2d.f90.s"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/n286654/Desktop/Incompact3d/decomp2d/decomp_2d.f90 -o CMakeFiles/decomp2d.dir/decomp_2d.f90.s

decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o: decomp2d/CMakeFiles/decomp2d.dir/flags.make
decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o: ../decomp2d/glassman.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/n286654/Desktop/Incompact3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/n286654/Desktop/Incompact3d/decomp2d/glassman.f90 -o CMakeFiles/decomp2d.dir/glassman.f90.o

decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/decomp2d.dir/glassman.f90.i"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/n286654/Desktop/Incompact3d/decomp2d/glassman.f90 > CMakeFiles/decomp2d.dir/glassman.f90.i

decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/decomp2d.dir/glassman.f90.s"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/n286654/Desktop/Incompact3d/decomp2d/glassman.f90 -o CMakeFiles/decomp2d.dir/glassman.f90.s

decomp2d/CMakeFiles/decomp2d.dir/io.f90.o: decomp2d/CMakeFiles/decomp2d.dir/flags.make
decomp2d/CMakeFiles/decomp2d.dir/io.f90.o: ../decomp2d/io.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/n286654/Desktop/Incompact3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object decomp2d/CMakeFiles/decomp2d.dir/io.f90.o"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/n286654/Desktop/Incompact3d/decomp2d/io.f90 -o CMakeFiles/decomp2d.dir/io.f90.o

decomp2d/CMakeFiles/decomp2d.dir/io.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/decomp2d.dir/io.f90.i"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/n286654/Desktop/Incompact3d/decomp2d/io.f90 > CMakeFiles/decomp2d.dir/io.f90.i

decomp2d/CMakeFiles/decomp2d.dir/io.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/decomp2d.dir/io.f90.s"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/n286654/Desktop/Incompact3d/decomp2d/io.f90 -o CMakeFiles/decomp2d.dir/io.f90.s

decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o: decomp2d/CMakeFiles/decomp2d.dir/flags.make
decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o: ../decomp2d/fft_generic.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/n286654/Desktop/Incompact3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/n286654/Desktop/Incompact3d/decomp2d/fft_generic.f90 -o CMakeFiles/decomp2d.dir/fft_generic.f90.o

decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/decomp2d.dir/fft_generic.f90.i"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/n286654/Desktop/Incompact3d/decomp2d/fft_generic.f90 > CMakeFiles/decomp2d.dir/fft_generic.f90.i

decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/decomp2d.dir/fft_generic.f90.s"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && /usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/n286654/Desktop/Incompact3d/decomp2d/fft_generic.f90 -o CMakeFiles/decomp2d.dir/fft_generic.f90.s

# Object files for target decomp2d
decomp2d_OBJECTS = \
"CMakeFiles/decomp2d.dir/decomp_2d.f90.o" \
"CMakeFiles/decomp2d.dir/glassman.f90.o" \
"CMakeFiles/decomp2d.dir/io.f90.o" \
"CMakeFiles/decomp2d.dir/fft_generic.f90.o"

# External object files for target decomp2d
decomp2d_EXTERNAL_OBJECTS =

lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/decomp_2d.f90.o
lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/glassman.f90.o
lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/io.f90.o
lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/fft_generic.f90.o
lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/build.make
lib/libdecomp2d.a: decomp2d/CMakeFiles/decomp2d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/n286654/Desktop/Incompact3d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking Fortran static library ../lib/libdecomp2d.a"
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && $(CMAKE_COMMAND) -P CMakeFiles/decomp2d.dir/cmake_clean_target.cmake
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/decomp2d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
decomp2d/CMakeFiles/decomp2d.dir/build: lib/libdecomp2d.a
.PHONY : decomp2d/CMakeFiles/decomp2d.dir/build

decomp2d/CMakeFiles/decomp2d.dir/clean:
	cd /home/n286654/Desktop/Incompact3d/build/decomp2d && $(CMAKE_COMMAND) -P CMakeFiles/decomp2d.dir/cmake_clean.cmake
.PHONY : decomp2d/CMakeFiles/decomp2d.dir/clean

decomp2d/CMakeFiles/decomp2d.dir/depend:
	cd /home/n286654/Desktop/Incompact3d/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/n286654/Desktop/Incompact3d /home/n286654/Desktop/Incompact3d/decomp2d /home/n286654/Desktop/Incompact3d/build /home/n286654/Desktop/Incompact3d/build/decomp2d /home/n286654/Desktop/Incompact3d/build/decomp2d/CMakeFiles/decomp2d.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : decomp2d/CMakeFiles/decomp2d.dir/depend


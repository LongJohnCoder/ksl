# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/roger/Projects/ksl

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/roger/Projects/ksl/src

# Include any dependencies generated for this target.
include test/CMakeFiles/check_vector.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/check_vector.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/check_vector.dir/flags.make

test/CMakeFiles/check_vector.dir/check_vector.c.o: test/CMakeFiles/check_vector.dir/flags.make
test/CMakeFiles/check_vector.dir/check_vector.c.o: ../test/check_vector.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/roger/Projects/ksl/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object test/CMakeFiles/check_vector.dir/check_vector.c.o"
	cd /home/roger/Projects/ksl/src/test && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/check_vector.dir/check_vector.c.o   -c /home/roger/Projects/ksl/test/check_vector.c

test/CMakeFiles/check_vector.dir/check_vector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/check_vector.dir/check_vector.c.i"
	cd /home/roger/Projects/ksl/src/test && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/roger/Projects/ksl/test/check_vector.c > CMakeFiles/check_vector.dir/check_vector.c.i

test/CMakeFiles/check_vector.dir/check_vector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/check_vector.dir/check_vector.c.s"
	cd /home/roger/Projects/ksl/src/test && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/roger/Projects/ksl/test/check_vector.c -o CMakeFiles/check_vector.dir/check_vector.c.s

test/CMakeFiles/check_vector.dir/check_vector.c.o.requires:

.PHONY : test/CMakeFiles/check_vector.dir/check_vector.c.o.requires

test/CMakeFiles/check_vector.dir/check_vector.c.o.provides: test/CMakeFiles/check_vector.dir/check_vector.c.o.requires
	$(MAKE) -f test/CMakeFiles/check_vector.dir/build.make test/CMakeFiles/check_vector.dir/check_vector.c.o.provides.build
.PHONY : test/CMakeFiles/check_vector.dir/check_vector.c.o.provides

test/CMakeFiles/check_vector.dir/check_vector.c.o.provides.build: test/CMakeFiles/check_vector.dir/check_vector.c.o


# Object files for target check_vector
check_vector_OBJECTS = \
"CMakeFiles/check_vector.dir/check_vector.c.o"

# External object files for target check_vector
check_vector_EXTERNAL_OBJECTS =

test/check_vector: test/CMakeFiles/check_vector.dir/check_vector.c.o
test/check_vector: test/CMakeFiles/check_vector.dir/build.make
test/check_vector: libksl.a
test/check_vector: /usr/lib/libm.so
test/check_vector: test/CMakeFiles/check_vector.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/roger/Projects/ksl/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable check_vector"
	cd /home/roger/Projects/ksl/src/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/check_vector.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/check_vector.dir/build: test/check_vector

.PHONY : test/CMakeFiles/check_vector.dir/build

test/CMakeFiles/check_vector.dir/requires: test/CMakeFiles/check_vector.dir/check_vector.c.o.requires

.PHONY : test/CMakeFiles/check_vector.dir/requires

test/CMakeFiles/check_vector.dir/clean:
	cd /home/roger/Projects/ksl/src/test && $(CMAKE_COMMAND) -P CMakeFiles/check_vector.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/check_vector.dir/clean

test/CMakeFiles/check_vector.dir/depend:
	cd /home/roger/Projects/ksl/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/roger/Projects/ksl /home/roger/Projects/ksl/test /home/roger/Projects/ksl/src /home/roger/Projects/ksl/src/test /home/roger/Projects/ksl/src/test/CMakeFiles/check_vector.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/check_vector.dir/depend


# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.16

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\Software\CLion 2019.3.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\Software\CLion 2019.3.3\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Project\DSA\as2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Project\DSA\as2\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/as2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/as2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/as2.dir/flags.make

CMakeFiles/as2.dir/main.cpp.obj: CMakeFiles/as2.dir/flags.make
CMakeFiles/as2.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Project\DSA\as2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/as2.dir/main.cpp.obj"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\as2.dir\main.cpp.obj -c D:\Project\DSA\as2\main.cpp

CMakeFiles/as2.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/as2.dir/main.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Project\DSA\as2\main.cpp > CMakeFiles\as2.dir\main.cpp.i

CMakeFiles/as2.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/as2.dir/main.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Project\DSA\as2\main.cpp -o CMakeFiles\as2.dir\main.cpp.s

# Object files for target as2
as2_OBJECTS = \
"CMakeFiles/as2.dir/main.cpp.obj"

# External object files for target as2
as2_EXTERNAL_OBJECTS =

as2.exe: CMakeFiles/as2.dir/main.cpp.obj
as2.exe: CMakeFiles/as2.dir/build.make
as2.exe: CMakeFiles/as2.dir/linklibs.rsp
as2.exe: CMakeFiles/as2.dir/objects1.rsp
as2.exe: CMakeFiles/as2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\Project\DSA\as2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable as2.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\as2.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/as2.dir/build: as2.exe

.PHONY : CMakeFiles/as2.dir/build

CMakeFiles/as2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\as2.dir\cmake_clean.cmake
.PHONY : CMakeFiles/as2.dir/clean

CMakeFiles/as2.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Project\DSA\as2 D:\Project\DSA\as2 D:\Project\DSA\as2\cmake-build-debug D:\Project\DSA\as2\cmake-build-debug D:\Project\DSA\as2\cmake-build-debug\CMakeFiles\as2.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/as2.dir/depend

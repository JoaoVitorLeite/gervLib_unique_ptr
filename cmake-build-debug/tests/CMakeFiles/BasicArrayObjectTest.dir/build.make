# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /snap/clion/248/bin/cmake/linux/x64/bin/cmake

# The command to remove a file.
RM = /snap/clion/248/bin/cmake/linux/x64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/joaovictor/Code/gervLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/joaovictor/Code/gervLib/cmake-build-debug

# Include any dependencies generated for this target.
include tests/CMakeFiles/BasicArrayObjectTest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/BasicArrayObjectTest.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/BasicArrayObjectTest.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/BasicArrayObjectTest.dir/flags.make

tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o: tests/CMakeFiles/BasicArrayObjectTest.dir/flags.make
tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o: /home/joaovictor/Code/gervLib/tests/BasicArrayObjectTest.cpp
tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o: tests/CMakeFiles/BasicArrayObjectTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/joaovictor/Code/gervLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o"
	cd /home/joaovictor/Code/gervLib/cmake-build-debug/tests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o -MF CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o.d -o CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o -c /home/joaovictor/Code/gervLib/tests/BasicArrayObjectTest.cpp

tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.i"
	cd /home/joaovictor/Code/gervLib/cmake-build-debug/tests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/joaovictor/Code/gervLib/tests/BasicArrayObjectTest.cpp > CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.i

tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.s"
	cd /home/joaovictor/Code/gervLib/cmake-build-debug/tests && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/joaovictor/Code/gervLib/tests/BasicArrayObjectTest.cpp -o CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.s

# Object files for target BasicArrayObjectTest
BasicArrayObjectTest_OBJECTS = \
"CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o"

# External object files for target BasicArrayObjectTest
BasicArrayObjectTest_EXTERNAL_OBJECTS =

tests/BasicArrayObjectTest: tests/CMakeFiles/BasicArrayObjectTest.dir/BasicArrayObjectTest.cpp.o
tests/BasicArrayObjectTest: tests/CMakeFiles/BasicArrayObjectTest.dir/build.make
tests/BasicArrayObjectTest: tests/CMakeFiles/BasicArrayObjectTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/joaovictor/Code/gervLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable BasicArrayObjectTest"
	cd /home/joaovictor/Code/gervLib/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BasicArrayObjectTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/BasicArrayObjectTest.dir/build: tests/BasicArrayObjectTest
.PHONY : tests/CMakeFiles/BasicArrayObjectTest.dir/build

tests/CMakeFiles/BasicArrayObjectTest.dir/clean:
	cd /home/joaovictor/Code/gervLib/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/BasicArrayObjectTest.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/BasicArrayObjectTest.dir/clean

tests/CMakeFiles/BasicArrayObjectTest.dir/depend:
	cd /home/joaovictor/Code/gervLib/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/joaovictor/Code/gervLib /home/joaovictor/Code/gervLib/tests /home/joaovictor/Code/gervLib/cmake-build-debug /home/joaovictor/Code/gervLib/cmake-build-debug/tests /home/joaovictor/Code/gervLib/cmake-build-debug/tests/CMakeFiles/BasicArrayObjectTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/BasicArrayObjectTest.dir/depend


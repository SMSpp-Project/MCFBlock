# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MCFBlock                                       #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
project(MCFBlock VERSION 0.3.1 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
include(CTest)

# This will set an user-settable cache variable.
# The value of the variable is propagated to consuming libraries via the
# *Config.cmake file.
set(MCFBlock_SOLVER "cplex" CACHE STRING "Solver used by MCFBlock.")
set_property(CACHE MCFBlock_SOLVER PROPERTY STRINGS relax cplex mfsmx sptree)

if (MCFBlock_SOLVER STREQUAL "relax")
    message(STATUS "MCFBlock: using RelaxIV solver")
    add_compile_definitions(HAVE_RELAX)
elseif (MCFBlock_SOLVER STREQUAL "cplex")
    message(STATUS "MCFBlock: using CPLEX solver")
    add_compile_definitions(HAVE_CPLEX)
elseif (MCFBlock_SOLVER STREQUAL "mfsmx")
    message(STATUS "MCFBlock: using Simplex solver")
    add_compile_definitions(HAVE_MFSMX)
elseif (MCFBlock_SOLVER STREQUAL "sptree")
    message(STATUS "MCFBlock: using Shortest Path Tree solver")
    add_compile_definitions(HAVE_SPTRE)
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)
find_package(MCFClass REQUIRED)

if (MCFBlock_SOLVER STREQUAL "cplex" AND NOT MCFClass_USE_CPLEX)
    message(FATAL_ERROR "MCFBlock solver is \"CPLEX\" but MCFCplex was not built")
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/MCFBlockConfig.h.in MCFBlockConfig.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(MCFBlock)
target_compile_features(MCFBlock PUBLIC cxx_std_17)

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
target_sources(MCFBlock PRIVATE
               src/MCFSolver.cpp
               src/MCFBlock.cpp)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path..
target_include_directories(
        MCFBlock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(MCFBlock PUBLIC
                      SMS++::SMSpp
                      MCFClass::MCFClass)

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(SMS++::MCFBlock ALIAS MCFBlock)

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS MCFBlock
        EXPORT MCFBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/MCFBlockConfig.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT MCFBlock-targets
        FILE MCFBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# ----- Export from the build tree ------------------------------------------ #
# The following commands are needed so other projects can find the library
# in a local build directory. They use CMake's User Package Registry
# (located in ~/.cmake/registry) to do so.
# See: https://cmake.org/cmake/help/latest/manual/cmake-packages.7.html
# They don't affect the installation.

# CMAKE_EXPORT_PACKAGE_REGISTRY is OFF by default.
# It may be ON if the project is included by the SMS++ umbrella project.
if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    # Generate the target import file, that allows other
    # CMake projects to import the target.
    export(EXPORT MCFBlock-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockTargets.cmake
           NAMESPACE SMS++::)

    # Export the package in user's package registry
    export(PACKAGE MCFBlock)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

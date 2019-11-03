# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MCF Block                                      #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8...3.15)
# TODO: set proper version
project(MCFBlock VERSION 1.0.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    enable_testing()
endif ()

# MCF Solver
set(MCF_SOLVER "relax" CACHE STRING "MCF Solver")
set_property(CACHE MCF_SOLVER PROPERTY STRINGS relax cplex mfsmx sptree)

if (MCF_SOLVER STREQUAL "relax")
    message(STATUS "Using RelaxIV solver")
    add_compile_definitions(HAVE_RELAX)
elseif (MCF_SOLVER STREQUAL "cplex")
    message(STATUS "Using CPLEX solver")
    add_compile_definitions(HAVE_CPLEX)
elseif (MCF_SOLVER STREQUAL "mfsmx")
    message(STATUS "Using Simplex solver")
    add_compile_definitions(HAVE_MFSMX)
elseif (MCF_SOLVER STREQUAL "sptree")
    message(STATUS "Using Shortest Path Tree solver")
    add_compile_definitions(HAVE_SPTRE)
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(SMSpp REQUIRED)
find_package(MCFClass REQUIRED)

# ----- Library ------------------------------------------------------------- #
add_library(MCFBlock)
target_compile_features(MCFBlock PUBLIC cxx_std_17)

target_sources(MCFBlock PRIVATE
               src/MCFSolver.cpp
               src/MCFBlock.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        MCFBlock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(MCFBlock PUBLIC
                      SMS++::SMSpp
                      MCFClass::MCFClass)

# Alias for testing
add_library(SMS++::MCFBlock ALIAS MCFBlock)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS MCFBlock
        EXPORT MCFBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT MCFBlock-targets
        FILE MCFBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# Install package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MCFBlock)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)

    export(EXPORT MCFBlock-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockTargets.cmake
           NAMESPACE SMS++::)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE MCFBlock)
endif ()

# ----- Tools --------------------------------------------------------------- #
add_subdirectory(tools)

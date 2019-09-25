# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MCF Block                                      #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8)
# ----- Project description and directories --------------------------------- #
# TODO: set proper version
project(MCFBlock VERSION 1.0.0)

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Requirements -------------------------------------------------------- #
# find_package() is used to find libraries and automatically configure them.
# It requires either a package configuration module (***Config.cmake) or a
# find module (Find***.cmake), usually provided by libraries.

# If we are configuring from the umbrella project, we assume that the valid
# settings are there.
# If we are building this project alone we use its settings but we check
# also for settings in its parent directory.
# TODO: Check if this behaviour is correct
if (NOT SMS++_CONFIG_ALL)
    include(CMakeCustom.txt OPTIONAL)
    include(../CMakeCustom.txt OPTIONAL)
endif ()

# This is to use the find modules we provide
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(SMSpp REQUIRED)
find_package(MCFClass REQUIRED)
find_package(netCDF REQUIRED)

# ----- Library ------------------------------------------------------------- #
add_library(MCFBlock STATIC)
target_compile_features(MCFBlock PUBLIC cxx_std_17)

target_sources(MCFBlock PRIVATE
               ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
               ${PROJECT_SOURCE_DIR}/MCFBlock.cpp
               )

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(MCFBlock PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>)

target_link_libraries(MCFBlock PUBLIC
                      SMS++::SMSpp
                      MCFClass::MCFClass
                      NetCDF::NetCDFCxx)

# Alias for testing
add_library(SMS++::MCFBlock ALIAS MCFBlock)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

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
        DESTINATION ${INSTALL_CONFIGDIR})

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
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install package version and configuration files
# TODO: Check if we should install also find modules
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        # ${CMAKE_CURRENT_LIST_DIR}/cmake/FindnetCDF.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

export(EXPORT MCFBlock-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockTargets.cmake
       NAMESPACE SMS++::)

# Export also the find modules for the local config file to use it
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindnetCDF.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindnetCDF.cmake
               COPYONLY)

# FIXME: I think that exporting this is a bad idea
#configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeCustom.txt
#               ${CMAKE_CURRENT_BINARY_DIR}/CMakeCustom.txt
#               COPYONLY)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MCFBlock)

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(test)
add_subdirectory(tools)

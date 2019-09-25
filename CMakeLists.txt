# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MCF Block                                      #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.8)
# ----- Project description and directories -----------------------------------
# TODO: set proper version
project(MCFBlock VERSION 1.0.0)

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Options for choosing the MCF Solver -----------------------------------
# TODO: Do it in a more fashionable way
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
#add_compile_definitions(HAVE_RELAX)
add_definitions(-DHAVE_RELAX)

# ----- Requirements ----------------------------------------------------------
# This file contains HINTS supported by the find modules
include(CMakeConfig.txt OPTIONAL)

find_package(SMSpp    REQUIRED)
find_package(MCFClass REQUIRED)
find_package(netCDF   REQUIRED)

# ----- Library ---------------------------------------------------------------
add_library(MCFBlock STATIC)
target_compile_features(MCFBlock PUBLIC cxx_std_17)

target_sources(MCFBlock PRIVATE
    ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
    ${PROJECT_SOURCE_DIR}/MCFBlock.cpp
    )

target_include_directories(MCFBlock PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>)

target_link_libraries(MCFBlock PUBLIC
                      SMS++::SMSpp
                      MCFClass::MCFClass
                      NetCDF::NetCDFCxx)

# Alias for testing
add_library(SMS++::MCFBlock ALIAS MCFBlock)

# ----- Install instructions --------------------------------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Install library
install(TARGETS MCFBlock
        EXPORT  MCFBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install export file *Targets.cmake
install(EXPORT MCFBlock-targets
        FILE   MCFBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${INSTALL_CONFIGDIR})

# Create *ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Create *Config.cmake
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR})

# Install the Config and ConfigVersion
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR})

# ----- Export from the build tree --------------------------------------------
export(EXPORT MCFBlock-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockTargets.cmake
       NAMESPACE SMS++::)

# Export the find modules and the config file so other projects will use them
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindnetCDF.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/FindnetCDF.cmake
               COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/CMakeConfig.txt
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeConfig.txt
               COPYONLY)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MCFBlock)

# ----- Subdirectories --------------------------------------------------------
add_subdirectory(test)
add_subdirectory(tools)

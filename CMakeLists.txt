# ----- Project description and directories ------------------------------------
cmake_minimum_required(VERSION 3.8)
project(MCFBlock
        VERSION 1.0.0
        DESCRIPTION "MCF Block for SMS++")

set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_SOURCE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----- Options for choosing the MCF Solver ------------------------------------
# TODO: Do it in a more fashionable way
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(SMSpp REQUIRED)
find_package(MCFClass REQUIRED)

# ----- Library ----------------------------------------------------------------
set(MCFBlock_SOURCES
    ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
    ${PROJECT_SOURCE_DIR}/MCFBlock.cpp
    )

# TODO: Find out why STATIC doesn't work
add_library(MCFBlock OBJECT ${MCFBlock_SOURCES})

# Alias for testing
add_library(SMS++::MCFBlock ALIAS MCFBlock)

target_compile_features(MCFBlock PUBLIC cxx_std_17)
target_include_directories(MCFBlock
                           PUBLIC
                           $<INSTALL_INTERFACE:include/SMS++>
                           $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>
                           )

target_link_libraries(MCFBlock PUBLIC
                      SMS++::SMSpp
                      MCFClass::MCFClass
                      Boost::boost
                      Eigen3::Eigen)
# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
add_subdirectory(tools)

# ----- Install instructions ---------------------------------------------------
include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SMS++)

# Install library
install(TARGETS MCFBlock
        EXPORT MCFBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

# Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Export target in *Targets.cmake
install(EXPORT MCFBlock-targets
        FILE MCFBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

# Create *ConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

# Create *Config.cmake
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/MCFBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Install the Config and ConfigVersion
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT MCFBlock-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/MCFBlockTargets.cmake
       NAMESPACE SMS++::)

# Register package in user's package registry (~/.cmake/registry)
export(PACKAGE MCFBlock)

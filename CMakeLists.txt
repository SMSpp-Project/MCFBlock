cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MCFBlock DESCRIPTION "MCF Block for SMS++")

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

# ----- Options for choosing the MCF Solver ------------------------------------
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# ----- Includes ---------------------------------------------------------------
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${PROJECT_INCLUDE_DIR})

# TODO: Write and use a FindSMS++.cmake
include_directories(../SMS++/include)

# TODO: Write and use Find*.cmake files
set(MCFClass_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../MCFClass)
include_directories(${MCFClass_INCLUDE_DIR}/MCFClass
                    ${MCFClass_INCLUDE_DIR}/MCFClone
                    ${MCFClass_INCLUDE_DIR}/MCFCplex
                    ${MCFClass_INCLUDE_DIR}/MCFSimplex
                    ${MCFClass_INCLUDE_DIR}/OPTUtils
                    ${MCFClass_INCLUDE_DIR}/RelaxIV
                    )

# ----- SOURCE FILES -----------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
    ${PROJECT_SOURCE_DIR}/MCFBlock.cpp
    )

# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
add_subdirectory(tools)

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(MCFBlock DESCRIPTION "MCF Block for SMS++")

# Options for choosing the MCF Solver
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

set(MCFBlock_SOURCES
    ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
    ${PROJECT_SOURCE_DIR}/MCFBlock.cpp)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# TODO: Write and use Find*.cmake files
set(MCFClass_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../MCFClass)
include_directories(${MCFClass_INCLUDE_DIR}/MCFClass
                    ${MCFClass_INCLUDE_DIR}/MCFClone
                    ${MCFClass_INCLUDE_DIR}/MCFCplex
                    ${MCFClass_INCLUDE_DIR}/MCFSimplex
                    ${MCFClass_INCLUDE_DIR}/OPTUtils
                    ${MCFClass_INCLUDE_DIR}/RelaxIV
                    )

include_directories(../SMS++/include)

include_directories(${PROJECT_INCLUDE_DIR})
# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${MCFBlock_SOURCES})

add_subdirectory(test)
add_subdirectory(tools)

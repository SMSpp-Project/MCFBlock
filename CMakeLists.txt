cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_VERBOSE_MAKEFILE OFF)

# ----- Project description and directories ------------------------------------
project(MCFBlock DESCRIPTION "MCF Block for SMS++")

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

# ----- Options for choosing the MCF Solver ------------------------------------
# TODO: Do it in a more fashionable way
#add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)
add_compile_definitions(HAVE_RELAX)

# ----- Requirements -----------------------------------------------------------
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# ----- Custom includes --------------------------------------------------------
# TODO: When libraries will be installable, fix this
set(SMS++_INCLUDE_DIR ../SMS++/include)
set(MCFClass_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../MCFClass)
set(MCFClass_INCLUDE_DIRS
    ${MCFClass_INCLUDE_DIR}/MCFClass
    ${MCFClass_INCLUDE_DIR}/MCFClone
    ${MCFClass_INCLUDE_DIR}/MCFCplex
    ${MCFClass_INCLUDE_DIR}/MCFSimplex
    ${MCFClass_INCLUDE_DIR}/OPTUtils
    ${MCFClass_INCLUDE_DIR}/RelaxIV
    )

# ----- Library ----------------------------------------------------------------
set(${PROJECT_NAME}_SOURCES
    ${PROJECT_SOURCE_DIR}/MCFSolver.cpp
    ${PROJECT_SOURCE_DIR}/MCFBlock.cpp
    )

# TODO: Find out why STATIC doesn't work
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${SMS++_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${MCFClass_INCLUDE_DIRS})

# ----- Subdirectories ---------------------------------------------------------
add_subdirectory(test)
add_subdirectory(tools)

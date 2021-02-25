# --------------------------------------------------------------------------- #
#    CMake file for unit tests                                                #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else ()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- MCFBlock_test ------------------------------------------------------- #
add_executable(MCFBlock_test test.cpp)
target_link_libraries(MCFBlock_test PRIVATE SMS++::MCFBlock)

# This version replicates the batch file in CMake.
# This method was preferred because it allows to keep track
# of the single tests via the CTest application and also
# avoids copying all the auxiliary files.
set(what "data/N3-5-5-2-2.nc4")
foreach (mode 0 1 2 12 13 14 20 29 30)
    foreach (seed RANGE 0 19)
        add_test(NAME "MCFBlock_test/[${what},${seed},${mode}]"
                 COMMAND MCFBlock_test ${what} ${seed} ${mode}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endforeach ()
endforeach ()

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#    CMake file for MCFBlock/test                                             #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (MSVC)
        add_link_options("/WHOLEARCHIVE")
        add_link_options("/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- Policies ------------------------------------------------------------ #
# Introduced in 3.19 to handle arbitrary characters in test names
# See: https://cmake.org/cmake/help/latest/policy/CMP0110.html
if (POLICY CMP0110)
    cmake_policy(SET CMP0110 NEW)
endif ()

# ----- MCFBlock_test ------------------------------------------------------- #
add_executable(MCFBlock_test test.cpp)
target_link_libraries(MCFBlock_test PRIVATE SMS++::MCFBlock)

# This version replicates the batch file in CMake.
# This method was preferred because it allows to keep track
# of the single tests via the CTest application and also
# avoids copying all the auxiliary files.
#set(what "data/N3-0-0-0-0.nc4")
#foreach (mode 0 1 2 12 13 14 20 29 30)
#    foreach (seed RANGE 0 19)
#        add_test(NAME "MCFBlock_test/[${what},${seed},${mode}]"
#                 COMMAND MCFBlock_test ${what} ${seed} ${mode}
#                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#    endforeach ()
#endforeach ()

# --------------------------------------------------------------------------- #

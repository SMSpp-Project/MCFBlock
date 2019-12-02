# --------------------------------------------------------------------------- #
#    CMake file for unit tests                                                #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from statically linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.

if (NOT BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- Legacy tests -------------------------------------------------------- #
add_executable(MCFBlock_test test.cpp)
target_link_libraries(MCFBlock_test PRIVATE SMS++::MCFBlock)

#add_test(
#        NAME MCFBlock_test
#        COMMAND MCFBlock_test
#)

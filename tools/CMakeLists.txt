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

# ----- dmx2nc4 ------------------------------------------------------------- #
add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_libraries(dmx2nc4 SMS++::MCFBlock)

# ----- Install instructions ------------------------------------------------ #
install(TARGETS dmx2nc4
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

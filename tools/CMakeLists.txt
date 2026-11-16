# --------------------------------------------------------------------------- #
#    CMake file for MCFBlock/tools                                            #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- dmx2nc4 ------------------------------------------------------------- #
add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_libraries(dmx2nc4 PRIVATE SMS++::MCFBlock)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)
install(TARGETS dmx2nc4
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --------------------------------------------------------------------------- #

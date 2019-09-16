add_executable(MCFBlock_test test.cpp)
target_link_libraries(MCFBlock_test
        SMS++
        MCFBlock
        MCFClass
        netcdf
        netcdf-cxx4)

add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_libraries(dmx2nc4
        SMS++
        MCFBlock
        MCFClass
        netcdf netcdf-cxx4)

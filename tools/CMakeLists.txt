find_package(netCDF REQUIRED)
find_package(MCFClass REQUIRED)
find_package(MCFBlock REQUIRED)

add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_directories(dmx2nc4 PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(dmx2nc4
                      SMS++::SMSpp
                      SMS++::MCFBlock
                      MCFClass::MCFClass
                      netcdf
                      netcdf-cxx4
                      )

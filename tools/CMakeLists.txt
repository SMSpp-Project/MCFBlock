find_package(netCDF REQUIRED)

add_executable(dmx2nc4 dmx2nc4.cpp)
target_include_directories(dmx2nc4 PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(dmx2nc4 PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(dmx2nc4
                      SMS++
                      MCFBlock
                      MCFClass
                      netcdf
                      netcdf-cxx4
                      )

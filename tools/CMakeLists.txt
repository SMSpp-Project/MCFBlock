find_package(netCDF REQUIRED)
include_directories(${netCDF_INCLUDE_DIR})
link_directories(${netCDF_LIB_DIR})

add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_libraries(dmx2nc4
                      SMS++
                      MCFBlock
                      MCFClass
                      netcdf
                      netcdf-cxx4)

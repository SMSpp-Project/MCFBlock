find_package(netCDF REQUIRED)

add_executable(MCFBlock_test test.cpp)
target_include_directories(MCFBlock_test PRIVATE ${netCDF_INCLUDE_DIR})
target_link_directories(MCFBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MCFBlock_test
                      SMS++
                      MCFBlock
                      MCFClass
                      netcdf
                      netcdf-cxx4
                      )

find_package(netCDF REQUIRED)
find_package(MCFClass REQUIRED)
find_package(MCFBlock REQUIRED)

add_executable(MCFBlock_test test.cpp)
target_link_directories(MCFBlock_test PRIVATE ${netCDF_LIB_DIR})
target_link_libraries(MCFBlock_test
                      SMS++::SMSpp
                      SMS++::MCFBlock
                      MCFClass::MCFClass
                      netcdf
                      netcdf-cxx4
                      )

add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)

add_executable(MCFBlock_test test.cpp)

target_link_libraries(MCFBlock_test
        SMS++
        MCFBlock
        MCFClass
#        MILPSolver
#        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libcplex.a
#        /Applications/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic/libilocplex.a
        netcdf
        netcdf-cxx4)

add_executable(dmx2nc4 dmx2nc4.cpp)
target_link_libraries(dmx2nc4
        SMS++
        MCFBlock
        MCFClass
        netcdf netcdf-cxx4)

add_compile_definitions(HAVE_MFSMX)
#add_compile_definitions(HAVE_CPLEX)

add_library(MCFBlock OBJECT MCFSolver.cpp MCFBlock.cpp)

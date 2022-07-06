whatis([[chrest-debug]])

help([[Loads latest debug version of Ablate and all required libraries.]])

family("chrest")

project_dir = "/projects/academic/chrest"
lib_dir = pathJoin(project_dir, "/lib")

setenv("LIB_DIR", lib_dir)
setenv("PROJECT_DIR", project_dir)

-- CCR modules
load("chrest/ccr")

--unload("intel-mpi/2020.2", "intel/20.2")

--load("openmpi/4.0.4")

-- CHREST modules
load("petsc-chrest/debug")
load("ablate/debug")
load("petscXdmfGenerator/v0.0.8-debug")

-- TChem modules
load("yaml-cpp")
load("gtest")
load("openblas/debug")
load("kokkos/debug")
load("tines/debug")
load("tchem/debug")

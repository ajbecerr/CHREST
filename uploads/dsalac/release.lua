whatis([[chrest]])

help([[Loads latest release version of Ablate and all required libraries.]])

family("chrest")

project_dir = "/projects/academic/chrest"
lib_dir = pathJoin(project_dir, "/lib")

setenv("LIB_DIR", lib_dir)
setenv("PROJECT_DIR", project_dir)

-- CCR modules
load("chrest/ccr")

-- CHREST modules
load("petsc-chrest/release")
load("ablate/release")
load("petscXdmfGenerator/v0.0.8-release")

-- TChem modules
load("yaml-cpp")
load("gtest")
load("openblas/release")
load("kokkos/release")
load("tines/release")
load("tchem/release")

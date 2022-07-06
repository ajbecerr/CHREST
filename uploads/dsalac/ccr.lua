whatis([[chrest_ccr]])

help([[Loads required CCR modules.]])

family("chrest_ccr")

-- CCR modules
load("intel/20.2")
load("intel-mpi/2020.2")
--load("intel/19.5", "intel-mpi/2019.5");
--load('hdf5/1.12.0-mpi');
load("gcc/11.2.0")
load("cmake/3.22.3")
load("valgrind/3.14.0")
load("gdb/7.8")


-- Some enviroment variables aren't set.
setenv("I_MPI_PMI_LIBRARY", "/usr/lib64/libpmi.so")
setenv("TEST_MPI_COMMAND", "srun")

mpi_bin = pathJoin(os.getenv("I_MPI_ROOT"), "intel64/bin")

setenv("CC", pathJoin(mpi_bin, "mpicc"))
setenv("CXX", pathJoin(mpi_bin, "mpicxx"))
setenv("FC", pathJoin(mpi_bin, "mpif90"))

--setenv("CC", "mpicc")
--setenv("CXX", "mpicxx")
--setenv("FC", "mpif90")

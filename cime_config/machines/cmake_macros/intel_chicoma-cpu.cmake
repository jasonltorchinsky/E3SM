set(PIO_FILESYSTEM_HINTS "lustre")
string(APPEND CONFIG_ARGS " --host=cray")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -qmkl")

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "icc")
set(SCXX "icpc")
set(SFC "ifc")

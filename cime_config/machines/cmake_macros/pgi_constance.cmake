string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CPPDEFS " -DLINUX")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -traceback  -fpe0 -check  all -check noarg_temp_created -ftrapuv")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -C -Mchkfpstk -Mchkstk -Mdalign  -Mdepchk -Miomutex -Mrecursive -Meh_frame")
set(PIO_FILESYSTEM_HINTS "lustre")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpmi -L$ENV{MPI_LIB} -lmpich")

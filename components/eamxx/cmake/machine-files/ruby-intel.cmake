include(${CMAKE_CURRENT_LIST_DIR}/quartz.cmake)
set(CMAKE_CXX_FLAGS "-w -cxxlib=/usr/tce/packages/gcc/gcc-8.3.1/rh" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/ -mkl" CACHE STRING "" FORCE)
set(PYTHON_EXECUTABLE "/usr/tce/packages/python/python-3.8.2/bin/python3" CACHE STRING "" FORCE)
set(RUN_ML_CORRECTION_TEST TRUE CACHE BOOL "")
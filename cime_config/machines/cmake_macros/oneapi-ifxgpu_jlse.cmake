
set(CXX_LINKER "CXX")
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nf-config --flibs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0} -Wl,-rpath -Wl,$ENV{NETCDF_PATH}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nc-config --libs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0}")
string(APPEND SLIBS " -fiopenmp -fopenmp-targets=spir64")
set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")
set(USE_SYCL "TRUE")
string(APPEND SYCL_FLAGS " -\-intel -fsycl -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\"") # for pvc node only
#string(APPEND SYCL_FLAGS " -\-intel -fsycl")
string(APPEND CXX_LDFLAGS " -Wl,-\-defsym,main=MAIN_\_ -lifcore -\-intel -fsycl -lsycl -Xsycl-target-backend \"-device 12.60.7\"") #for pvc node only

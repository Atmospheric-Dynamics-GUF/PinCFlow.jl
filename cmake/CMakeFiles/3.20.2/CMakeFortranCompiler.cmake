set(CMAKE_Fortran_COMPILER "/sw/spack-levante/intel-oneapi-mpi-2021.5.0-efceba/mpi/2021.5.0/bin/mpiifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "2021.5.0.20211109")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/sw/spack-levante/intel-oneapi-mpi-2021.5.0-efceba/mpi/2021.5.0/include;/sw/spack-levante/intel-oneapi-compilers-2022.0.1-an2cbq/compiler/2022.0.1/linux/compiler/include/intel64;/sw/spack-levante/intel-oneapi-compilers-2022.0.1-an2cbq/compiler/2022.0.1/linux/compiler/include/icc;/sw/spack-levante/intel-oneapi-compilers-2022.0.1-an2cbq/compiler/2022.0.1/linux/compiler/include;/usr/local/include;/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0/include;/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0/include-fixed;/sw/spack-levante/gcc-11.2.0-bcn7mb/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;dl;rt;pthread;ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/sw/spack-levante/intel-oneapi-mpi-2021.5.0-efceba/mpi/2021.5.0/lib/release;/sw/spack-levante/intel-oneapi-mpi-2021.5.0-efceba/mpi/2021.5.0/lib;/sw/spack-levante/intel-oneapi-compilers-2022.0.1-an2cbq/compiler/2022.0.1/linux/compiler/lib/intel64_lin;/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0;/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64;/lib64;/usr/lib64;/sw/spack-levante/gcc-11.2.0-bcn7mb/lib;/lib;/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

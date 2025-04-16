# CMakeLists of PinCFlow

# Specify minimum cmake version.
cmake_minimum_required(VERSION 2.8)

# Set project name and language.
project(pinc Fortran)

# Set directories for object files and executable.
set(CMAKE_Fortran_MODULE_DIRECTORY "../build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")

# Set compiler.
set(CMAKE_Fortran_COMPILER mpif90)
# set(CMAKE_Fortran_COMPILER mpiifort)

# Set debug option.
#option(DEBUGFLAG "" OFF)
 option(DEBUGFLAG "" ON)

# Set compiler options.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(DEBUGFLAG)
    add_compile_options(-Og -g -fcheck=all -Wall -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra -Wno-unused-variable -fallow-argument-mismatch)
  else()
    add_compile_options(-O3 -fdefault-real-8 -fbacktrace -funroll-loops -fallow-argument-mismatch -w)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if(DEBUGFLAG)
    add_compile_options(-O0 -g -fpe0 -check all -real-size 64 -traceback -unroll=4 -ip)
  else()
    add_compile_options(-O3 -real-size 64 -traceback -unroll=4 -ip)
  endif()
endif()

# Set netCDF flags.
execute_process(COMMAND bash "-c" "echo `nf-config --fflags --flibs` -Wl,-rpath,`nf-config --prefix`/lib" OUTPUT_VARIABLE NCFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Add netCDF flags.
set(CMAKE_Fortran_COMPILE_OBJECT "${CMAKE_Fortran_COMPILE_OBJECT} ${NCFLAGS}")
set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} ${NCFLAGS}")

# Add executable.
add_executable(pinc
    ../src/atmosphere.f90
    ../src/bicgstab_tools.f90
    ../src/boundary.f90
    ../src/finish.f90
    ../src/fluxes.f90
    ../src/ice.f90
    ../src/init.f90
    ../src/mpi.f90
    ../src/muscl.f90
    ../src/output_netcdf.f90
    ../src/pinc.f90
    ../src/poisson.f90
    ../src/sizeof.f90
    ../src/timeScheme.f90
    ../src/tracer.f90
    ../src/types.f90
    ../src/update.f90
    ../src/wkb.f90)

# Print information.
message("[INFO] Compiler ID is " ${CMAKE_Fortran_COMPILER_ID})
message("[INFO] DEBUGFLAG is " ${DEBUGFLAG})

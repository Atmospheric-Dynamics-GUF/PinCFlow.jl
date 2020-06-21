cmake_minimum_required(VERSION 2.8)
project(pinc Fortran)

# set Fortran compiler
# set(CMAKE_Fortran_COMPILER /usr/bin/mpif90)
set(MPI_Fortran_COMPILER mpif90)

# set language to fortran
enable_language(Fortran)

# set directories for output binary independent of IDE definitions
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_Fortran_MODULE_DIRECTORY build)

# set compiler options
add_compile_options(-O0 -g -fcheck=all -Wall -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra -Wno-unused-variable)

# find MPI module
find_package( MPI )
include_directories( ${MPI_Fortran_INCLUDE_PATH} )

# compile executable
add_executable(pinc
        src/types.f90
        src/atmosphere.f90
        src/mpi.f90
        src/timeScheme.f90
        src/algebra.f90
        src/init.f90
        src/debug.f90
        src/muscl.f90
        src/wkb.f90
        src/xweno.f90
        src/fluxes.f90
        src/boundary.f90
        src/hypre_tools.f90
        src/poisson.f90
        src/update.f90
        src/output.f90
        src/finish.f90
        src/pinc.f90
        src/ice.f90
        )

# link MPI and hypre from system libraries
target_link_libraries( pinc ${MPI_Fortran_LIBRARIES} HYPRE)

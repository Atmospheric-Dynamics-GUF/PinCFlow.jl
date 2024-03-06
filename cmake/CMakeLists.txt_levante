# CMakeLists of PincFlowMSGWAM

# Specify minimum cmake version.
cmake_minimum_required(VERSION 2.8)

# Set compiler.
# set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_Fortran_COMPILER mpiifort)

# Set project name.
project(pinc Fortran)

# Set language.
enable_language(Fortran)

# Set debugq option.
option(DEBUGFLAG "set to ON to enable debug compiler options" OFF)

# Set directories for output binary independent of IDE definitions.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
set(CMAKE_Fortran_MODULE_DIRECTORY "../build")

# Set compiler options.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (DEBUGFLAG)
        add_compile_options(-O0 -g -fcheck=all -Wall -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra -Wno-unused-variable)
    else()
        add_compile_options(-O3 -fdefault-real-8 -fbacktrace -funroll-loops)
        # add_compile_options(-O3 -fdefault-real-8 -fbacktrace -funroll-loops -fallow-argument-mismatch)
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (DEBUGFLAG)
        add_compile_options(-O0 -g -march=core-avx2 -fp-model source -fpe0 -check all -real-size 64 -traceback -unroll=4 -fnoinline -no-ip )
    else()
        add_compile_options(-O3 -march=core-avx2 -fp-model source -real-size 64 -traceback -unroll=4 -ip)
    endif()
endif()

# Find MPI module.
find_package(MPI)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Set link to HYPRE library.
# set(HYPRE /usr/lib/x86_64-linux-gnu/libHYPRE-2.22.1.so)
# set(HYPRE /home/atmodynamics/jochum/libraries/spack/opt/spack/linux-scientific7-haswell/gcc-4.8.5/hypre-2.26.0-fgfu4ncdxl7ullqryad7jtdzgzqfacjw/lib/libHYPRE.so)

# Compile executable.
add_executable(pinc
        ../src/types.f90
        ../src/sizeof.f90
        ../src/mpi.f90
        ../src/timeScheme.f90
        ../src/algebra.f90
        ../src/atmosphere.f90
        ../src/init.f90
        ../src/debug.f90
        ../src/muscl.f90
        ../src/wkb.f90
        ../src/xweno.f90
        ../src/fluxes.f90
        ../src/boundary.f90
        ../src/bicgstab_tools.f90
        ../src/poisson.f90
        ../src/update.f90
        ../src/output.f90
        ../src/finish.f90
        ../src/ice.f90
        ../src/tracer.f90
        ../src/pinc.f90
        )

# Link MPI and Hypre from system libraries.
target_link_libraries(pinc ${MPI_Fortran_LIBRARIES})# ${HYPRE})

# Inform about debug options.
message("[INFO] Compiler ID is " ${CMAKE_Fortran_COMPILER_ID})
message("[INFO] DEBUGFLAG is " ${DEBUGFLAG})
add_custom_command(TARGET pinc POST_BUILD COMMAND echo "[INFO] DEBUGFLAG is " ${DEBUGFLAG})


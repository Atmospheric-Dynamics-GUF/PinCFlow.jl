# Installing NetCDF-Fortran and Required Libraries

## Introduction

> NetCDF (Network Common Data Form) is a set of software libraries and machine-independent data formats. It supports the creation, access, and sharing of array-oriented scientific data. NetCDF is widely used in scientific applications and is a community standard for data sharing. The Unidata Program Center supports netCDF programming interfaces for C, C++, Java, Fortran, and others.

On the Levante cluster, `netcdf-fortran` is already installed. You simply need to load the module before compiling and running PinCFlow.

On the Goethe Cluster, however, you will need to install several libraries manually: zlib, HDF5, netCDF-C, and netCDF-Fortran. The latter three libraries must be installed with parallel support enabled. This guide walks you through installing these libraries step by step.

Refer to the official documentation for more details: [NetCDF C Library](https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html) and [NetCDF Fortran Library](https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html).

## Step-by-Step Guide

### 1. Download the Required Libraries

Use `wget` to download the source files for the following libraries:

- [zlib](https://www.zlib.net/)
- [HDF5](https://portal.hdfgroup.org/downloads/)
- [netCDF-C](https://downloads.unidata.ucar.edu/netcdf/)
- [netCDF-Fortran](https://downloads.unidata.ucar.edu/netcdf/)

For example, you can download `zlib` with:

```bash
wget https://www.zlib.net/zlib-1.3.1.tar.gz # Link from September 13, 2024
```

### 2. Extract the Downloaded Files

Extract the downloaded `.tar.gz` or `.zip` files in your preferred directory:

```bash
tar -xvzf source-file.tar.gz
```

or

```bash
unzip source-file.zip
```

### 3. Create an Installation Directory

Create a directory where you want to install the libraries. For example:

```bash
mkdir /home/atmodynamics/name/my-libraries
```

### 4. Set Up Environment Variables for Compiling

Set the necessary environment variables for compilation. Make sure to adjust these paths and flags according to your specific setup:

```bash
export F77=`which mpif77`
export FC=`which mpif90`
export CC=`which mpicc`
export CXX=`which mpicxx`
export CFLAGS=-DgFortran
export LIBDIR=/home/atmodynamics/name/my-libraries
```

*These settings are for the GNU compiler. If you're using the Intel compiler, change the compiler flags accordingly (e.g., use ifort, icc, etc.).*

### 5. Install zlib

Move to the `zlib` source folder and configure the installation:

```bash
cd zlib-<version>
./configure --prefix=${LIBDIR}
make check
make install
```

### 6. Install HDF5 with Parallel Support

Move to the `hdf5` source folder and install it with parallel support enabled:

```bash
cd hdf5-<version>
./configure --prefix=${LIBDIR} --with-zlib=${LIBDIR} --with-default-plugindir=${LIBDIR}/plugin --enable-parallel
make check
make install
```

### 7. Install NetCDF-C with Parallel Support

Move to the `netcdf-c` source folder and install it with parallel support:

```bash
cd netcdf-c-<version>
LDFLAGS=-L${LIBDIR}/lib CPPFLAGS=-I${LIBDIR}/include ./configure --prefix=${LIBDIR} --with-plugin-dir=${LIBDIR}/plugin
make check
make install
```

### 8. Install NetCDF-Fortran with Parallel Support

Move to the `netcdf-fortran` source folder and install it:

```bash
cd netcdf-fortran-<version>
LDFLAGS=-L${LIBDIR}/lib CPPFLAGS=-I${LIBDIR}/include LD_LIBRARY_PATH=${LIBDIR}/lib HDF5_PLUGIN_PATH=${LIBDIR}/plugin LIBS="-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl" ./configure --prefix=${LIBDIR} --disable-shared
make check
make install
```

### 9. Update Environment Variables

To ensure your system can locate the libraries, update your `.bashrc` file with the following environment variables. This will ensure they are available in future sessions:

```bash
export NETCDF=/home/atmodynamics/name/my-libraries
export PATH=$NETCDF/bin:$PATH
export NETCDF_INCDIR=$NETCDF/include
export NETCDF_LIBDIR=$NETCDF/lib
export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:$NETCDF/lib:$LD_LIBRARY_PATH
export PATH NETCDF
```

After updating your `.bashrc`, run:

```bash
source ~/.bashrc
```

### 10. Compile and Link Programs with NetCDF Fortran Libraries

Now that everything is installed, you can compile and link your Fortran program with the NetCDF libraries. Use the following command to compile:

```bash
mpif90 -o my_prog my_prog.f90 `nf-config --fflags --flibs` -Wl,-rpath,`nf-config --prefix`/lib
```

This command uses `nf-config` to automatically provide the necessary compiler flags for linking with the NetCDF Fortran libraries.

## Additional Notes

- Ensure that you have the necessary permissions to install libraries in the chosen directories.
- Use the `make check` command to verify that everything was installed correctly.
- If you encounter any issues, check the log files generated during the `make` steps for troubleshooting information.
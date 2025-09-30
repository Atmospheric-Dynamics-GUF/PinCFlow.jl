# PinCFlow.jl

## Introduction

PinCFlow is an atmospheric flow solver that was designed for performing idealized simulations. It integrates the Boussinesq, pseudo-incompressible and compressible equations in a conservative flux form ([Klein, 2009](https://doi.org/10.1007/s00162-009-0104-y); [Rieper et al., 2013](https://doi.org/10.1175/mwr-d-12-00026.1)), using a a semi-implicit method that combines explicit and implicit time-stepping schemes ([Benacchio & Klein, 2019](https://doi.org/10.1175/MWR-D-19-0073.1); [Schmid et al., 2021](https://doi.org/10.1175/MWR-D-21-0126.1)). Spatially, the equations are discretized with a finite-volume method, such that all quantities are represented by averages over grid cells and fluxes are computed on the respective cell interfaces. The grid is staggered so that the velocity components are defined at the same points as the corresponding fluxes of scalar quantities. PinCFlow operates in a vertically stretched terrain-following coordinate system based on [Gal-Chen and Somerville (1975a)](https://doi.org/10.1016/0021-9991(75)90037-6), [Gal-Chen and Somerville (1975b)](https://doi.org/10.1016/0021-9991(75)90054-6) and [Clark (1977)](https://doi.org/10.1016/0021-9991(77)90057-2).

The Lagrangian WKB model MSGWaM is interactively coupled to PinCFlow, so that unresolved gravity waves may be parameterized in a manner that accounts for transience and horizontal propagation. The resolved fields are updated with tendencies computed by MSGWaM at the beginning of every time step. A description of the theory behind PinCFlow-MSGWaM can be found in [Achatz et al. (2017)](https://doi.org/10.1002/qj.2926) and [Achatz et al. (2023)](https://doi.org/10.1063/5.0165180). For a numerical perspective and more information on the development, see [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381), [Boeloeni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1), [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1), [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1) and [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1).

## User guide

### Installation

To install PinCFlow, first make sure you have installed [Julia](https://docs.julialang.org/en/v1/manual/installation/). You can then clone this repository with

```shell
git clone git@github.com:Atmospheric-Dynamics-GUF/PinCFlow.jl.git
```

and set up the project environment by running

```shell
julia --project -e 'using Pkg; Pkg.instantiate()'
```

in the root directory of your clone.

### Running the model

As a minimal example, the script

```julia
using PinCFlow

integrate(Namelists())
```

runs PinCFlow in its default configuration, if executed with

```shell
julia --project script.jl
```

in the root directory of the repository. This simulation will finish comparatively quickly and won't produce particularly interesting results, since PinCFlow simply initializes a $1 \times 1 \times 1 \, \mathrm{km^3}$ isothermal atmosphere at rest with $3 \times 3 \times 3$ grid points and integrates the governing equations over one hour. A more complex configuration can be set up by providing namelists with changed parameters. For instance, running the script

```julia
# examples/submit/periodic_hill.jl

using PinCFlow

@ivy if length(ARGS) == 0
    output_file = "./pincflow_output.h5"
elseif length(ARGS) == 1
    output_file = ARGS[1] * "/pincflow_output.h5"
else
    error("Too many arguments to the script!")
end

atmosphere = AtmosphereNamelist(; initial_wind = (1.0E+1, 0.0E+0, 0.0E+0))
domain = DomainNamelist(;
    ndx = 40,
    ndy = 1,
    ndz = 40,
    lx = 2.0E+4,
    ly = 2.0E+4,
    lz = 2.0E+4,
)
grid = GridNamelist(; mountainheight_dim = 1.0E+1, mountainwidth_dim = 1.0E+4)
output = OutputNamelist(; output_variables = (:w,), output_file = output_file)
sponge = SpongeNamelist(; spongelayer = true)

integrate(Namelists(; atmosphere, domain, grid, output, sponge))

```

yields a 2D simulation with an initial wind of $10 \, \mathrm{m \, s^{- 1}}$ that generates a mountain wave above a periodic hill. The vertical wind is written to the output file `pincflow_output.h5` in the directory specified by an additional argument to the script (or the current directory, if that argument is omitted). More involved examples are given in the "Examples" section of the documentation. A description of all namelists and their parameters is provided in the "Reference" section.

If you want to run PinCFlow in parallel, make sure you are using the correct backends for [MPI.jl](https://juliaparallel.org/MPI.jl/latest/) and [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/). By default, the two packages use JLL backends that have been automatically installed. If you want to keep this setting, you only need to make sure to use the correct MPI binary (specifically not that of a default MPI installation on your system). You can do so by running

```shell
mpiexec=$(julia --project -e 'using MPICH_jll; println(MPICH_jll.mpiexec_path)')
${mpiexec} -n ${tasks} julia --project script.jl
```

with `tasks` set to the number of MPI processes. Note that in `script.jl`, the parameters `npx`, `npy` and `npz` of the namelist `domain`, which represent the number of MPI processes in the three dimensions of physical space, need to be set such that their product is equal to `tasks`.

However, if you plan to run PinCFlow on a cluster, you may want to consider using a provided MPI installation as backend. In that case, the MPI preferences need to be updated accordingly and the HDF5 backend has to be set to a library that has been installed with parallel support, using the chosen MPI installation. This can be done by running

```shell
julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names = ["/path/to/mpi/library/"])'
julia --project -e 'using HDF5; HDF5.API.set_libraries!("/path/to/libhdf5.so", "/path/to/libhdf5_hl.so")'
```

with the paths set appropriately (more details can be found in the documentations of MPI.jl and HDF5.jl). Note that this configuration will be saved in `LocalPreferences.toml`, so that the new backends will be used by all future scripts run in the project. By running

```shell
julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
julia --project -e 'using HDF5; HDF5.API.set_libraries!()'
```

you can restore the default backends. Having configured MPI.jl and HDF5.jl to use installations on your system, you can run

```shell
mpiexec -n ${tasks} julia --project script.jl
```

with `mpiexec` being your chosen system binary. For users who would like to run PinCFlow on [Goethe](https://csc.uni-frankfurt.de/wiki/doku.php?id=public:usage:goethe) or [Levante](https://docs.dkrz.de/doc/levante/index.html), shell-script examples are provided in `examples/submit`.

### Visualizing the results

PinCFlow uses parallel HDF5 to write simulation data. By default, the path to the output file is `pincflow_output.h5` (from the directory in which the run script is executed). This may be changed by setting the parameter `output_file` of the namelist `output` accordingly. The dimensions of most output fields are (in order) $\widehat{x}$ (zonal axis), $\widehat{y}$ (meridional axis), $\widehat{z}$ (axis orthogonal to the vertical coordinate surfaces) and $t$ (time). Ray-volume property fields differ slightly in that they have an additional (spectral) dimension in front and a vertical dimension that includes the first ghost layer below the surface. To specify which fields are to be written, set the parameters `output_variables`, `save_ray_volumes` and `prepare_restart` of the namelist `output` accordingly (more details are given in the "Reference" section of the documentation).

For the visualization of simulation results, we recommend using [PythonPlot.jl](https://github.com/JuliaPy/PythonPlot.jl). A function that configures PythonPlot.jl to use a preset style, as well as one that facilitates the generation of symmetric contour plots, are exported by `PinCFlow`. The script

```julia
# examples/visualization/periodic_hill.jl

using HDF5
using PythonPlot
using LaTeXStrings
using PinCFlow

set_plot_style()

# Import the data.
@ivy if length(ARGS) == 0
    data = h5open("./pincflow_output.h5")
elseif length(ARGS) == 1
    data = h5open(ARGS[1] * "/pincflow_output.h5")
else
    error("Too many arguments to the script!")
end

# Set the grid.
x = data["x"][:] ./ 1000
z = data["z"][:, 1, :] ./ 1000
x = x .* ones(size(z))

# Get the vertical wind.
w = data["w"][:, 1, :, end]

# Close the file.
close(data)

# Create the plot.
(levels, colormap) = symmetric_contours(minimum(w), maximum(w))
contours = contourf(x, z, w; levels = levels, cmap = colormap)
xlabel(L"x\,\left[\mathrm{km}\right]")
ylabel(L"z\,\left[\mathrm{km}\right]")
colorbar(contours; label = L"w\,\left[\mathrm{m\,s^{-1}}\right]")
savefig("examples/results/periodic_hill.png")
clf()

```

is an example for how to visualize the vertical wind at the end of a simple mountain-wave simulation performed with the script introduced above. Once again, the directory which the output file has been saved to is given as an additional argument to the script. The resulting plot is displayed below.

![](examples/results/periodic_hill.png)

## List of publications

 1. Initial flow solver: [Rieper et al. (2013)](https://doi.org/10.1175/mwr-d-12-00026.1)

 1. Initial gravity-wave scheme: [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381)

 1. Gravity-wave breaking scheme: [Boeloeni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1)

 1. Gravity-wave theory: [Achatz et al. (2017)](https://doi.org/10.1002/qj.2926)

 1. Coupling of the flow solver and gravity-wave scheme: [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1)

 1. Horizontal propagation and direct approach in the gravity-wave scheme: [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1)

 1. Semi-implicit time scheme: [Schmid et al. (2021)](https://doi.org/10.1175/MWR-D-21-0126.1)

 1. Extended gravity-wave theory: [Achatz et al. (2023)](https://doi.org/10.1063/5.0165180)

 1. Terrain-following coordinates & orographic source: [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1)

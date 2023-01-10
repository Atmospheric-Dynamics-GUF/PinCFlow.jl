# PincFlow

## Introduction

PincFlow solves the pseudo-incompressible equations in a conservative flux form, using either a third-order accurate Runge-Kutta method [(Rieper et al., 2013)](https://doi.org/10.1175/mwr-d-12-00026.1) or a semi-implicit method that combines explicit and implicit time-stepping schemes [(Schmid et al., 2021)](https://doi.org/10.1175/MWR-D-21-0126.1). The equations are discretized with a finite-volume method, such that all quantities are represented by spatial averages over grid cells and fluxes are computed on the respective cell interfaces. The grid is staggered so that the velocity components are defined at the same points as the corresponding fluxes of scalar quantities.

In an adiabatic configuration, the equations can also be solved above uneven ground, in which case PincFlow uses a terrain-following grid based on [Gal-Chen and Somerville (1975)](https://doi.org/10.1016/0021-9991(75)90037-6). The topography can be set by specifying the corresponding namelist parameters. Topographic data can also be read from an appropriate input file (the topographic output of the model can be used as a reference).

The Lagrangian WKB model MS-GWaM is coupled interactively to PincFlow, such that unresolved gravity waves may be parameterized. The resolved fields are then updated according to the tendencies computed by the ray tracer at every Runge-Kutta substep. A description of MS-GWaM can be found in [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381), [Bölöni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1) and [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1).

## Code organization

In addition to the source code (`src`), the following resources are provided.

* A `Makefile` that can be used to compile the code with either `mpif90` or `mpiifort`, as well as a corresponding `CMakeLists.txt` file (`cmake`)

* Namelist files, run scripts, visualization tools and sample plots for a set of canonical test cases (`tests`)

* A code formatter that can be used to unify spacing, indentation and linebreaks of all Fortan files in a given directory (`tools`)

## Workflow

The code is shared in a GitLab repository. Any contributions to the code should adhere to the following workflow.

1. A development branch should be created from the master branch or synchronized accordingly. **No direct changes are to be made on the master branch.**

2. If the contributor does not already have a connected local repository, the remote repository is to be cloned.

3. A local branch should be created from the remote development branch. The indended changes to the code are to be made on this branch.

4. The changes may be pushed to the remote development branch, e.g. to let other contributors review them.

5. The remote master branch must be updated regularly. For this purpose, remote development branches are merged with it. **Before any merge can be done, the corresponding development branch must be stable and all canonical test cases must run successfully on it. The merge itself is to be checked by at least one additional contributor with sufficient experience.**

Whenever the master branch is updated, all contributors should synchronize their development branches correspondingly.
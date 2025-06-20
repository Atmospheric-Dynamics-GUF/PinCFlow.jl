# PinCFlow

## Introduction

PinCFlow integrates the pseudo-incompressible equations in a conservative flux form ([Klein, 2009](https://doi.org/10.1007/s00162-009-0104-y); [Rieper et al., 2013](https://doi.org/10.1175/mwr-d-12-00026.1)), using a a semi-implicit method that combines explicit and implicit time-stepping schemes ([Benacchio & Klein, 2019](https://doi.org/10.1175/MWR-D-19-0073.1); [Schmid et al., 2021](https://doi.org/10.1175/MWR-D-21-0126.1)). The equations are discretized with a finite-volume method, such that all quantities are represented by spatial averages over grid cells and fluxes are computed on the respective cell interfaces. The grid is staggered so that the velocity components are defined at the same points as the corresponding fluxes of scalar quantities. PinCFlow operates in a vertically stretched terrain-following coordinate system based on [Gal-Chen and Somerville (1975a)](https://doi.org/10.1016/0021-9991(75)90037-6), [Gal-Chen and Somerville (1975b)](https://doi.org/10.1016/0021-9991(75)90054-6) and [Clark (1977)](https://doi.org/10.1016/0021-9991(77)90057-2).

The Lagrangian WKB model MS-GWaM is interactively coupled to PinCFlow, so that unresolved gravity waves may be parameterized in a manner that accounts for transience and horizontal propagation. The resolved fields are updated with tendencies computed by MS-GWaM at the beginning of every time step. A description of PinCFlow-MS-GWaM can be found in [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1), [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1) and [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1).

## Workflow

The code is shared in a GitLab repository. Any contributions to the code should adhere to the following workflow.

1. If you are new to the project, create a remote development branch for your contributions (name it such that others can identify it as your branch) and clone the repository.

1. Make your changes on your local development branch.

1. Pull recent changes made on the remote master branch into your local master branch and merge it into your local development branch, resolving merge conflicts if necessary.

1. **Ensure that the model is stable and that all canonical tests reproduce the sample results.**

1. Push your changes to your remote development branch.

1. Request to merge your remote development branch into the remote master branch.

## List of publications

1. [Rieper et al. (2013)](https://doi.org/10.1175/mwr-d-12-00026.1)
1. [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1)
1. [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1)
1. [Schmid et al. (2021)](https://doi.org/10.1175/MWR-D-21-0126.1)
1. [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1)

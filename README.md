# PinCFlow

## Introduction

PinCFlow integrates the pseudo-incompressible equations in a conservative flux form, using either a third-order accurate Runge-Kutta method [(Rieper et al., 2013)](https://doi.org/10.1175/mwr-d-12-00026.1) or a semi-implicit method that combines explicit and implicit time-stepping schemes [(Schmid et al., 2021)](https://doi.org/10.1175/MWR-D-21-0126.1). The equations are discretized with a finite-volume method, such that all quantities are represented by spatial averages over grid cells and fluxes are computed on the respective cell interfaces. The grid is staggered so that the velocity components are defined at the same points as the corresponding fluxes of scalar quantities. PinCFlow operates in a vertically stretched terrain-following coordinate system based on [Gal-Chen and Somerville (1975a)](https://doi.org/10.1016/0021-9991(75)90037-6), [Gal-Chen and Somerville (1975b)](https://doi.org/10.1016/0021-9991(75)90054-6) and [Clark (1977)](https://doi.org/10.1016/0021-9991(77)90057-2).

## Workflow

The code is shared in a GitLab repository. Any contributions to the code should adhere to the following workflow.

1. Create a fork from the upstream repository.

1. Clone the fork to create a local repository.

1. Make your changes on the development branch of the local repository.

1. Commit your changes in reproducible steps.

1. Pull the development branch of the upstream repository and merge it into the development branch of the local repository, resolving merge conflicts if necessary.

1. **Ensure that the model is stable and that all canonical tests reproduce the sample results.**

1. Push the development branch of the local repository to the fork.

1. Request to merge the development branch of the fork into the development branch of the upstream repository.

![](workflow.png "Workflow  ")
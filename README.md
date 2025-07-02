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

## Building and accessing the documentation

The code uses Documenter.jl to build the documentation. To build the documentation, run the following command in the root directory of the repository:

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The documentation will be generated in the `docs/build` directory. To view the documentation, open the `index.html` file in a web browser.

## Style guide

### Code

1. Put every module, composite type (including constructor methods) and function into a file on its own, with the file name matching that of the object. Create a folder for every module.

1. Only use Unicode in strings and docstrings.

### Documentation

1. Write a docstring for every module, type and method.

1. Always show the exact full signature of the object at the top of the documentation, within a Julia code block.

1. Include a single sentence describing what the function does or what the object represents after the signature block. Use the imperative form when documenting methods. If needed, provide more details in a second paragraph, after a blank line.

1. List all arguments with descriptions (but without types and default values) under an `# Arguments` header, with one `-` bullet for each argument.

1. If a method returns something other than `nothing`, list all returned objects with descriptions under a `# Returns` header, with one `-` bullet for each returned object. Include the types (in Julia syntax).

1. Provide links to related functions under a `See also` header, with one `-` bullet for each function.

1. Use single backticks to identify code and double backticks to identify equations. Use Unicode characters rather than LaTeX escape sequences.

1. Place the starting and ending `"""` characters on lines by themselves.

## List of publications

1. [Rieper et al. (2013)](https://doi.org/10.1175/mwr-d-12-00026.1)
1. [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1)
1. [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1)
1. [Schmid et al. (2021)](https://doi.org/10.1175/MWR-D-21-0126.1)
1. [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1)

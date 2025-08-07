# PinCFlow

## Introduction

PinCFlow integrates the pseudo-incompressible equations in a conservative flux form ([Klein, 2009](https://doi.org/10.1007/s00162-009-0104-y); [Rieper et al., 2013](https://doi.org/10.1175/mwr-d-12-00026.1)), using a a semi-implicit method that combines explicit and implicit time-stepping schemes ([Benacchio & Klein, 2019](https://doi.org/10.1175/MWR-D-19-0073.1); [Schmid et al., 2021](https://doi.org/10.1175/MWR-D-21-0126.1)). The equations are discretized with a finite-volume method, such that all quantities are represented by spatial averages over grid cells and fluxes are computed on the respective cell interfaces. The grid is staggered so that the velocity components are defined at the same points as the corresponding fluxes of scalar quantities. PinCFlow operates in a vertically stretched terrain-following coordinate system based on [Gal-Chen and Somerville (1975a)](https://doi.org/10.1016/0021-9991(75)90037-6), [Gal-Chen and Somerville (1975b)](https://doi.org/10.1016/0021-9991(75)90054-6) and [Clark (1977)](https://doi.org/10.1016/0021-9991(77)90057-2).

The Lagrangian WKB model MSGWaM is interactively coupled to PinCFlow, so that unresolved gravity waves may be parameterized in a manner that accounts for transience and horizontal propagation. The resolved fields are updated with tendencies computed by MSGWaM at the beginning of every time step. A description of PinCFlow-MSGWaM can be found in [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1), [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1) and [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1).

## User guide

## Developer guide

### Workflow

The code is shared in a GitLab repository. Any contributions to the code should adhere to the following workflow.

1. If you are new to the project, create a remote development branch for your contributions (name it such that others can identify it as your branch) and clone the repository.

1. Make your changes on your local development branch.

1. Pull recent changes made on the remote master branch into your local master branch and merge it into your local development branch, resolving merge conflicts if necessary.

1. **Ensure that the model is stable and that all canonical tests reproduce the sample results.**

1. Push your changes to your remote development branch.

1. Request to merge your remote development branch into the remote master branch.

### Writing code

* Put every module, composite type (including constructor methods) and function into a file on its own, with the file name matching that of the object. Create a folder for every module.

* Do not use Unicode.

* Use `CamelCase` for the names of modules and types. Use single captial letters for type parameters. For all other objects, use `snake_case` (in case the name only contains (preferrably whole) words, e.g. `vertical_wind`) and `squashedcase` (in case the name is mathematical, e.g. `what` for $\widehat{w}$).

* Use parametric composite types.

* Declare the types of method arguments.

* Use `@views` for expressions that create slices.

### Writing documentation

* Write a docstring for every module, function and type.

* Module docstrings:

    1. Include the exact full signature within a Julia code block, followed by a single descriptive (pseudo-)sentence and (if needed) a second paragraph with more details.

    1. List links to imported modules in a `# See also` section, with one `-` bullet for each.

* Function docstrings:

    1. For every method, include the exact full signature within a Julia code block, followed by a single, descriptive sentence in imperative form and (if needed) a second paragraph with more details.

    1. List all positional and optional arguments with descriptions (but without types and default values) in an `# Arguments` section, with one `-` bullet for each.

    1. List all keyword arguments with descriptions (but without types and default values) in a `# Keywords` section, with one `-` bullet for each.

    1. If the methods of a function return something other than `nothing`, list all returned objects with descriptions in a `# Returns` section, with one `-` bullet for each.

    1. List links to constructors/functions that are called in any of the explicitly defined constructor methods in a # See also section, with one `-` bullet for each.

* Type docstrings:

    1. Include the exact full signature within a Julia code block, followed by a single descriptive (pseudo-)sentence and (if needed) a second paragraph with more details.

    1. If the type is composite, include the exact full signature within a Julia code block, followed by a single, descriptive sentence in imperative form and (if needed) a second paragraph with more details, for each explicitly defined constructor method.

    1. If the type is composite, list all fields with their type restrictions and descriptions in a `# Fields` section, with one `-` bullet for each.

    1. If the type is composite, list all positional and optional arguments of the explicitly defined constructor methods with descriptions (but without types and default values) in an `# Arguments` section, with one `-` bullet for each.

    1. If the type is composite, list all keyword arguments of the explicitly defined constructor methods with descriptions (but without types and default values) in a `# Keywords` section, with one `-` bullet for each.

    1. If the type is composite, list links to constructors/functions that are called in any of the explicitly defined constructor methods in a # See also section, with one `-` bullet for each.

* Use single backticks to identify code and double backticks to identify equations. Use LaTeX escape sequences rather than Unicode characters.

* Place the starting and ending `"""` characters on lines by themselves.

### Building and accessing the documentation

The code uses Documenter.jl to build the documentation. To build the documentation, run the following command in the root directory of the repository:

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The documentation will be generated in the `docs/build` directory. To view the documentation, open the `index.html` file in a web browser.

## List of publications

1. [Rieper et al. (2013)](https://doi.org/10.1175/mwr-d-12-00026.1)
1. [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381)
1. [Boeloeni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1)
1. [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1)
1. [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1)
1. [Schmid et al. (2021)](https://doi.org/10.1175/MWR-D-21-0126.1)
1. [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1)

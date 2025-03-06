# PinCFlow.jl

In order to run the code, please install `julia` https://julialang.org/downloads/. We strongly recommend using `juliaup` which is also explained in this link. This code has been tested with version 1.10.8, so please use that version to ensure reproducibility. Once you have `julia`, clone this repository, enter the current directory and start `julia` in the current environment as

```shell
julia --project=.
```
This is equivalent to starting `julia` and then entering `import Pkg; Pkg.activate(".")`. Now install the required dependencies by running

```julia
julia> import Pkg; Pkg.instantiate()
```
This is only needed to be done once or after some packages have been updated. You can run an available test case as
```julia
julia> include("examples/elixir_mountainwave.jl")
```
This will also create a `mountainwave.png` file for visualization. You can change the parameters used in the simulation through the elixir file.
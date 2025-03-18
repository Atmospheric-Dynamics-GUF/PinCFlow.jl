if isdefined(@__MODULE__, :LanguageServer)
  include("../src/PinCFlow.jl")
end

using Test
using PinCFlow
using Trixi: CompressibleEulerEquations1D
using LinearAlgebra
using GZip

# include("test_poisson.jl")
# include("test_elixirs.jl")
# include("test_fluxes.jl")
# include("test_update.jl")
include("test_atmosphere.jl")
# include("test_sponge.jl")
# include("test_boundary.jl")

# include("test_namelist.jl")
include("test_init.jl")

"""
```julia
FoundationalTypes
```

Module that contains the composite types `Time`, `Constants`, `Domain`, `Grid`, `Atmosphere` and `Sponge`.

# Scopes

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)

# Exports

  - [`PinCFlow.Types.FoundationalTypes.Atmosphere`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Constants`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Domain`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Grid`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Sponge`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.Time`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_halos_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_halos_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_halos_of_field!`](@ref)
"""
module FoundationalTypes

using MPI
using ..NamelistTypes
using ...Macros

include("Time.jl")
include("Constants.jl")
include("Domain.jl")
include("Grid.jl")
include("Atmosphere.jl")
include("Sponge.jl")

include("compute_topography.jl")
include("compute_n2!.jl")

include("../../Boundaries/set_zonal_boundaries_of_field!.jl")
include("../../Boundaries/set_meridional_boundaries_of_field!.jl")
include("../../Boundaries/set_vertical_boundaries_of_field!.jl")
include("../../MPIOperations/set_zonal_halos_of_field!.jl")
include("../../MPIOperations/set_meridional_halos_of_field!.jl")
include("../../MPIOperations/set_vertical_halos_of_field!.jl")

export Time, Constants, Domain, Grid, Atmosphere, Sponge

export set_zonal_boundaries_of_field!,
    set_meridional_boundaries_of_field!,
    set_vertical_boundaries_of_field!,
    set_zonal_halos_of_field!,
    set_meridional_halos_of_field!,
    set_vertical_halos_of_field!

end

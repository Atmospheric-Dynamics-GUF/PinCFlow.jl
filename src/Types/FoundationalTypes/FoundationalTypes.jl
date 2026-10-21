"""
```julia
FoundationalTypes
```

Module that contains the composite types `Time`, `Constants`, `Domain`, `Grid`, `Atmosphere` and `Sponge`.

# See also

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types.NamelistTypes`](@ref)
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

include("compute_new_grid!.jl")
include("compute_n2!.jl")
include("compute_background_fields!.jl")

include("../../Boundaries/set_zonal_boundaries_of_field!.jl")
include("../../Boundaries/set_meridional_boundaries_of_field!.jl")
include("../../Boundaries/set_vertical_boundaries_of_field!.jl")
include("../../MPIOperations/set_zonal_halos_of_field!.jl")
include("../../MPIOperations/set_meridional_halos_of_field!.jl")
include("../../MPIOperations/set_vertical_halos_of_field!.jl")

export Time, Constants, Domain, Grid, Atmosphere, Sponge

export compute_new_grid!,
    compute_background_fields!,
    set_zonal_boundaries_of_field!,
    set_meridional_boundaries_of_field!,
    set_vertical_boundaries_of_field!,
    set_zonal_halos_of_field!,
    set_meridional_halos_of_field!,
    set_vertical_halos_of_field!

end

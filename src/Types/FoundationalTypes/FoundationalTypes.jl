"""
```julia
FoundationalTypes
```

Module that contains the composite types `Time`, `Constants`, `Domain`, `Grid`, `Atmosphere` and `Sponge`.

# See also

  - [`PinCFlow.Types.NamelistTypes`](@ref)
"""
module FoundationalTypes

using MPI
using ..NamelistTypes

include("Time.jl")
include("Constants.jl")
include("Domain.jl")
include("Grid.jl")
include("Atmosphere.jl")
include("Sponge.jl")

include("compute_topography.jl")

# This is not ideal - maybe it would be better to have copies of these
# functions...
include("../../Boundaries/set_zonal_boundaries_of_field!.jl")
include("../../Boundaries/set_meridional_boundaries_of_field!.jl")
include("../../Boundaries/set_vertical_boundaries_of_field!.jl")
include("../../MPIOperations/set_zonal_halos_of_field!.jl")
include("../../MPIOperations/set_meridional_halos_of_field!.jl")
include("../../MPIOperations/set_vertical_halos_of_field!.jl")

export Time, Constants, Domain, Grid, Atmosphere, Sponge

end

module Boundaries

using ..Types
using ..MPIOperations

struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end
struct BoundaryGWIntegrals end
struct BoundaryGWTendencies end

include("set_boundaries!.jl")
include("set_meridional_boundaries_of_field!.jl")
include("set_meridional_boundaries_of_reduced_field!.jl")
include("set_meridional_boundaries!.jl")
include("set_vertical_boundaries!.jl")
include("set_zonal_boundaries_of_field!.jl")
include("set_zonal_boundaries_of_reduced_field!.jl")
include("set_zonal_boundaries!.jl")

export BoundaryPredictands,
    BoundaryReconstructions,
    BoundaryFluxes,
    BoundaryGWIntegrals,
    BoundaryGWTendencies

export set_boundaries!,
    set_meridional_boundaries_of_field!,
    set_meridional_boundaries_of_reduced_field!,
    set_meridional_boundaries!,
    set_vertical_boundaries!,
    set_zonal_boundaries_of_field!,
    set_zonal_boundaries_of_reduced_field!,
    set_zonal_boundaries!

end

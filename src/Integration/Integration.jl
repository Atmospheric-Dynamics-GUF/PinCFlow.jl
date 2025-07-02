module Integration

using MPI
using Dates
using ..Types
using ..Boundaries
using ..Update
using ..PoissonSolver
using ..FluxCalculator
using ..Output
using ..MSGWaM

include("compute_time_step.jl")
include("integrate.jl")
include("modify_compressible_wind!.jl")
include("reset_predictands!.jl")
include("save_backups!.jl")
include("synchronize_compressible_atmosphere!.jl")
include("synchronize_density_fluctuations!.jl")
include("explicit_integration!.jl")
include("implicit_integration!.jl")
include("wkb_integration!.jl")

export integrate

end

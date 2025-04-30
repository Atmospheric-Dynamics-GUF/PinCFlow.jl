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
    include("reset_fluxes!.jl")
    include("reset_predictands!.jl")
    include("synchronize_compressible_atmosphere!.jl")
    include("synchronize_density_fluctuations!.jl")

    export integrate

end

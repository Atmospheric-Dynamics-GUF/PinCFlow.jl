"""
```julia
PinCFlow
```

Main module of PinCFlow.jl.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Integration`](@ref)
"""
module PinCFlow

using MPI

include("@ivy.jl")
include("reduce_exceptions.jl")
include("ensemble.jl")
include("plot_output.jl")
include("set_visualization_theme!.jl")
include("symmetric_contours.jl")

export @ivy
export reduce_exceptions,
    ensemble, plot_output, set_visualization_theme!, symmetric_contours

include("Types/Types.jl")
include("MPIOperations/MPIOperations.jl")
include("Boundaries/Boundaries.jl")
include("Update/Update.jl")
include("PoissonSolver/PoissonSolver.jl")
include("FluxCalculator/FluxCalculator.jl")
include("Output/Output.jl")
include("MSGWaM/MSGWaM.jl")
include("Integration/Integration.jl")
include("Examples/Examples.jl")

using PrecompileTools
using .Types
using .Integration
using .Examples

@setup_workload begin
    redirect_stdio(; stderr = devnull, stdout = devnull) do
        mktempdir() do directory
            x_size = 3
            y_size = 3
            z_size = 5

            keywords = (
                output_file = directory * "pincflow_output.h5",
                visualize = false,
            )

            @compile_workload begin
                cold_bubble(; x_size, z_size, keywords...)
                hot_bubble(; x_size, z_size, keywords...)
                mountain_wave(; x_size, y_size, z_size, keywords...)
                periodic_hill(; x_size, z_size, keywords...)
                vortex(; x_size, y_size, keywords...)
                wave_packet(; x_size, y_size, z_size, keywords...)
                wkb_mountain_wave(; x_size, y_size, z_size, keywords...)
                wkb_wave_packet(; x_size, y_size, z_size, keywords...)
                wp_3d(; x_size, y_size, z_size, keywords...)
            end
            return
        end
        return
    end
end

# Export namelists.
export DomainNamelist,
    OutputNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    TurbulenceNamelist,
    Namelists

# Export model-state constructor.
export State

# Export integration function.
export integrate

# Export example functions.
export cold_bubble,
    hot_bubble,
    mountain_wave,
    periodic_hill,
    vortex,
    wave_packet,
    wkb_mountain_wave,
    wkb_wave_packet,
    wp_3d

end

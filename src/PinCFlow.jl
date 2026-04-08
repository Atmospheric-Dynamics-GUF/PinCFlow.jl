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

include("@ivy.jl")
include("plot_output.jl")
include("set_visualization_theme!.jl")
include("symmetric_contours.jl")

export @ivy
export plot_output, set_visualization_theme!, symmetric_contours

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
            x_size = 5
            y_size = 5
            z_size = 5

            output_file = directory * "pincflow_output.h5"

            visualize = false

            @compile_workload begin
                cold_bubble(; x_size, z_size, output_file, visualize)
                hot_bubble(; x_size, z_size, output_file, visualize)
                mountain_wave(; x_size, y_size, z_size, output_file, visualize)
                periodic_hill(; x_size, z_size, output_file, visualize)
                vortex(; x_size, y_size, output_file, visualize)
                wave_packet(; x_size, y_size, z_size, output_file, visualize)
                wkb_mountain_wave(;
                    x_size,
                    y_size,
                    z_size,
                    output_file,
                    visualize,
                )
                wkb_wave_packet(;
                    x_size,
                    y_size,
                    z_size,
                    output_file,
                    visualize,
                )
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
    wkb_wave_packet

end

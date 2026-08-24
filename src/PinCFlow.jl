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

function plot_output end
function set_visualization_theme! end
function symmetric_contours end

export plot_output, set_visualization_theme!, symmetric_contours

include("Macros/Macros.jl")
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
    redirect_stdout(devnull) do
        mktempdir() do directory
            keywords = (
                output_file = directory * "/pincflow_output.h5",
                visualize = false,
                x_size = 3,
                y_size = 3,
                z_size = 5,
            )

            for name in names(Examples)
                example = getproperty(Examples, name)
                if example isa Function
                    @compile_workload begin
                        example(; keywords...)
                    end
                end
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
    wave_Boussinesq,
    wkb_wave_Boussinesq

end

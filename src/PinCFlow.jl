"""
```julia
PinCFlow
```

Main module of PinCFlow.

# See also

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Integration`](@ref)
"""
module PinCFlow

using PythonPlot

include("@ivy.jl")
include("set_plot_style.jl")
include("symmetric_contours.jl")

export @ivy
export set_plot_style, symmetric_contours

include("Types/Types.jl")
include("MPIOperations/MPIOperations.jl")
include("Boundaries/Boundaries.jl")
include("Update/Update.jl")
include("PoissonSolver/PoissonSolver.jl")
include("FluxCalculator/FluxCalculator.jl")
include("Output/Output.jl")
include("MSGWaM/MSGWaM.jl")
include("Integration/Integration.jl")

using .Types
using .Integration

# Export namelists.
export DomainNamelist,
    OutputNamelist,
    SettingNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    Namelists

# Export singletons needed in namelists.
export Boussinesq, PseudoIncompressible, Compressible
export MountainWave, WKBMountainWave
export MCVariant
export UniformBoussinesq, StratifiedBoussinesq, Isothermal
export ExponentialSponge, COSMOSponge, PolynomialSponge, SinusoidalSponge
export ConstantWaveAction, ConstantWaveEnergy
export Box, Shapiro
export SteadyState, SingleColumn, MultiColumn
export NoTracer, LinearTracer

# Export model-state constructor.
export State

# Export integration function.
export integrate

# Used for testing
function examples_dir()
    return pkgdir(PinCFlow, "examples/submit")
end

export examples_dir

end

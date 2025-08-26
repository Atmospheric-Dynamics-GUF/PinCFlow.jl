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
    WavePacketNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    IceNamelist,
    TurbulenceNamelist,
    Namelists

# Export singletons needed in namelists.
export Boussinesq, PseudoIncompressible, Compressible
export MountainWave, WKBMountainWave, WavePacket
export PeriodicBoundaries, SolidWallBoundaries
export MCVariant
export UniformBoussinesq, StratifiedBoussinesq, Isothermal
export ExponentialSponge, COSMOSponge, PolynomialSponge, SinusoidalSponge
export ConstantWaveAction, ConstantWaveEnergy
export Box, Shapiro
export SteadyState, SingleColumn, MultiColumn
export NoTracer, LinearTracer
export NoIce, IceOn
export NoTurbulence, TurbulenceOn

# Export model-state constructor.
export State

# Export integration function.
export integrate

end

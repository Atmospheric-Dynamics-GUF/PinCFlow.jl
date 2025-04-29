module PinCFlow

using Dates

include("Types/Types.jl")
include("MPIOperations/MPIOperations.jl")
include("Boundaries/Boundaries.jl")
include("Update/Update.jl")
include("PoissonSolver/PoissonSolver.jl")
include("FluxCalculator/FluxCalculator.jl")
include("Output/Output.jl")
include("MSGWaM/MSGWaM.jl")

using .Types
using .Boundaries
using .Update
using .PoissonSolver
using .FluxCalculator
using .Output
using .MSGWaM

include("integrate.jl")

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
    Namelists

# Export singletons needed in namelists.
export Boussinesq, PseudoIncompressible, Compressible
export MountainWave, WKBMountainWave
export PeriodicBoundaries, SolidWallBoundaries
export MCVariant
export UniformBoussinesq, StratifiedBoussinesq, Isothermal
export FPlane
export ExponentialSponge, COSMOSponge, PolynomialSponge, SinusoidalSponge
export ConstantWaveAction, ConstantWaveEnergy
export Box, Shapiro
export SteadyState, SingleColumn, MultiColumn
export Clip, Scale

# Export integration function.
export integrate, State

end

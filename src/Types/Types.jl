module Types

abstract type AbstractVariable end
abstract type AbstractIntegration end

struct Rho <: AbstractVariable end
struct RhoP <: AbstractVariable end
struct U <: AbstractVariable end
struct V <: AbstractVariable end
struct W <: AbstractVariable end
struct PiP <: AbstractVariable end
struct EXPL <: AbstractIntegration end
struct IMPL <: AbstractIntegration end

include("NamelistTypes/NamelistTypes.jl")
include("FoundationalTypes/FoundationalTypes.jl")
include("PoissonTypes/PoissonTypes.jl")
include("VariableTypes/VariableTypes.jl")
include("WKBTypes/WKBTypes.jl")

using .NamelistTypes
using .FoundationalTypes
using .PoissonTypes
using .VariableTypes
using .WKBTypes

include("State.jl")

export AbstractBackground,
    AbstractCoriolis,
    AbstractLimiter,
    AbstractVariable,
    AbstractModel,
    AbstractTestCase,
    AbstractBoundaries,
    AbstractSponge,
    AbstractMergeMode,
    AbstractLaunchAlgorithm,
    AbstractIntegration,
    AbstractWKBMode,
    AbstractWKBTestCase,
    AbstractWKBFilter

export Total,
    Horizontal,
    Isothermal,
    ConstantCoriolis,
    MCVariant,
    Rho,
    RhoP,
    U,
    V,
    W,
    EXPL, IMPL,
    PiP,
    PseudoIncompressible,
    MountainWave,
    WKBMountainWave,
    PeriodicBoundaries,
    SolidWallBoundaries,
    ExponentialSponge,
    COSMOSponge,
    PolynomialSponge,
    SinusoidalSponge,
    ConstantWaveAction,
    ConstantWaveEnergy,
    Clip,
    Scale,
    SteadyState,
    SingleColumn,
    MultiColumn,
    Box,
    Shapiro

export DomainNamelist,
    OutputNamelist,
    SettingNamelist,
    DiscretizationNamelist,
    PoissonNamelist,
    AtmosphereNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    Namelists,
    Time,
    Constants,
    Domain,
    Grid,
    Atmosphere,
    Sponge,
    Tensor,
    Operator,
    Preconditioner,
    BicGStab,
    Correction,
    Poisson,
    Predictands,
    Tendencies,
    Backups,
    Auxiliaries,
    Reconstructions,
    Fluxes,
    Variables,
    GWIntegrals,
    GWTendencies,
    Rays,
    Increments,
    SurfaceIndices,
    WKB,
    State

end

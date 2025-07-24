module Types

abstract type AbstractVariable end

struct Rho <: AbstractVariable end
struct RhoP <: AbstractVariable end
struct U <: AbstractVariable end
struct V <: AbstractVariable end
struct W <: AbstractVariable end
struct PiP <: AbstractVariable end
struct P <: AbstractVariable end
struct Theta <: AbstractVariable end
struct Explicit end
struct Implicit end

include("NamelistTypes/NamelistTypes.jl")
include("FoundationalTypes/FoundationalTypes.jl")
include("PoissonTypes/PoissonTypes.jl")
include("VariableTypes/VariableTypes.jl")
include("WKBTypes/WKBTypes.jl")
include("TracerTypes/TracerTypes.jl")
include("IceTypes/IceTypes.jl")
include("TurbulenceTypes/TurbulenceTypes.jl")

using .NamelistTypes
using .FoundationalTypes
using .PoissonTypes
using .VariableTypes
using .WKBTypes
using .TracerTypes
using .IceTypes
using .TurbulenceTypes

include("State.jl")

export AbstractBackground,
    AbstractCoriolisMode,
    AbstractLimiter,
    AbstractVariable,
    AbstractModel,
    AbstractTestCase,
    AbstractBoundaries,
    AbstractSponge,
    AbstractMergeMode,
    AbstractWKBMode,
    AbstractWKBTestCase,
    AbstractWKBFilter,
    AbstractTracer,
    AbstractIce,
    AbstractTurbulence

export Rho,
    RhoP,
    U,
    V,
    W,
    PiP,
    P,
    Theta,
    Explicit,
    Implicit,
    UniformBoussinesq,
    StratifiedBoussinesq,
    Isothermal,
    FPlane,
    MCVariant,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
    MountainWave,
    WKBMountainWave,
    WavePacket,
    WKBWavePacket,
    PeriodicBoundaries,
    SolidWallBoundaries,
    ExponentialSponge,
    COSMOSponge,
    PolynomialSponge,
    SinusoidalSponge,
    ConstantWaveAction,
    ConstantWaveEnergy,
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
    WavePacketNamelist,
    GridNamelist,
    SpongeNamelist,
    WKBNamelist,
    TracerNamelist,
    IceNamelist,
    TurbulenceNamelist,
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
    Tracer,
    Ice,
    Turbulence,
    State,
    NoTracer,
    LinearTracer,
    LikeDensity,
    TracerPredictands,
    TracerAuxiliaries,
    TracerTendencies,
    TracerReconstructions,
    TracerForcings,
    TracerGWImpact,
    TracerFluxes,
    IceOn,
    NoIce,
    IcePredictands,
    IceAuxiliaries,
    IceTendencies,
    IceReconstructions,
    IceFluxes,
    TurbulenceOn,
    NoTurbulence,
    TurbulencePredictands,
    TurbulenceAuxiliaries,
    TurbulenceTendencies,
    TurbulenceReconstructions,
    TurbulenceFluxes

end

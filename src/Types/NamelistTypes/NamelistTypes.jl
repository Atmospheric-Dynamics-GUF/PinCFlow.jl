module NamelistTypes

abstract type AbstractBackground end
abstract type AbstractCoriolisMode end
abstract type AbstractLimiter end
abstract type AbstractModel end
abstract type AbstractTestCase end
abstract type AbstractBoundaries end
abstract type AbstractSponge end
abstract type AbstractMergeMode end
abstract type AbstractLaunchAlgorithm end
abstract type AbstractWKBMode end
abstract type AbstractWKBTestCase <: AbstractTestCase end
abstract type AbstractWKBFilter end

struct UniformBoussinesq <: AbstractBackground end
struct StratifiedBoussinesq <: AbstractBackground end
struct Isothermal <: AbstractBackground end
struct FPlane <: AbstractCoriolisMode end
struct MCVariant <: AbstractLimiter end
struct Boussinesq <: AbstractModel end
struct PseudoIncompressible <: AbstractModel end
struct Compressible <: AbstractModel end
struct MountainWave <: AbstractTestCase end
struct WKBMountainWave <: AbstractWKBTestCase end
struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end
struct ExponentialSponge <: AbstractSponge end
struct COSMOSponge <: AbstractSponge end
struct PolynomialSponge <: AbstractSponge end
struct SinusoidalSponge <: AbstractSponge end
struct ConstantWaveAction <: AbstractMergeMode end
struct ConstantWaveEnergy <: AbstractMergeMode end
struct Clip <: AbstractLaunchAlgorithm end
struct Scale <: AbstractLaunchAlgorithm end
struct SteadyState <: AbstractWKBMode end
struct SingleColumn <: AbstractWKBMode end
struct MultiColumn <: AbstractWKBMode end
struct Box <: AbstractWKBFilter end
struct Shapiro <: AbstractWKBFilter end

include("DomainNamelist.jl")
include("OutputNamelist.jl")
include("SettingNamelist.jl")
include("DiscretizationNamelist.jl")
include("PoissonNamelist.jl")
include("AtmosphereNamelist.jl")
include("GridNamelist.jl")
include("SpongeNamelist.jl")
include("WKBNamelist.jl")
include("Namelists.jl")

export AbstractBackground,
    AbstractCoriolisMode,
    AbstractLimiter,
    AbstractModel,
    AbstractTestCase,
    AbstractBoundaries,
    AbstractSponge,
    AbstractMergeMode,
    AbstractLaunchAlgorithm,
    AbstractWKBMode,
    AbstractWKBTestCase,
    AbstractWKBFilter

export UniformBoussinesq,
    StratifiedBoussinesq,
    Isothermal,
    FPlane,
    MCVariant,
    Boussinesq,
    PseudoIncompressible,
    Compressible,
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
    Namelists

end

module NamelistTypes

abstract type AbstractBackground end
abstract type AbstractCoriolisMode end
abstract type AbstractLimiter end
abstract type AbstractModel end
abstract type AbstractTestCase end
abstract type AbstractBoundaries end
abstract type AbstractSponge end
abstract type AbstractMergeMode end
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
struct WavePacket <: AbstractTestCase end
struct WKBMountainWave <: AbstractWKBTestCase end
struct WKBWavePacket <: AbstractWKBTestCase end
struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end
struct ExponentialSponge <: AbstractSponge end
struct COSMOSponge <: AbstractSponge end
struct PolynomialSponge <: AbstractSponge end
struct SinusoidalSponge <: AbstractSponge end
struct ConstantWaveAction <: AbstractMergeMode end
struct ConstantWaveEnergy <: AbstractMergeMode end
struct SteadyState <: AbstractWKBMode end
struct SingleColumn <: AbstractWKBMode end
struct MultiColumn <: AbstractWKBMode end
struct Box <: AbstractWKBFilter end
struct Shapiro <: AbstractWKBFilter end

abstract type AbstractTracer end
struct NoTracer <: AbstractTracer end
struct LinearTracer <: AbstractTracer end

abstract type AbstractIce end
struct NoIce <: AbstractIce end
struct IceOn <: AbstractIce end

abstract type AbstractTurbulence end
struct NoTurbulence <: AbstractTurbulence end
struct TurbulenceOn <: AbstractTurbulence end

include("DomainNamelist.jl")
include("OutputNamelist.jl")
include("SettingNamelist.jl")
include("DiscretizationNamelist.jl")
include("PoissonNamelist.jl")
include("AtmosphereNamelist.jl")
include("WavePacketNamelist.jl")
include("GridNamelist.jl")
include("SpongeNamelist.jl")
include("WKBNamelist.jl")
include("TracerNamelist.jl")
include("IceNamelist.jl")
include("TurbulenceNamelist.jl")
include("Namelists.jl")

export AbstractBackground,
    AbstractCoriolisMode,
    AbstractLimiter,
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
    Shapiro,
    NoTracer,
    LinearTracer,
    NoIce,
    IceOn,
    NoTurbulence,
    TurbulenceOn

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

end

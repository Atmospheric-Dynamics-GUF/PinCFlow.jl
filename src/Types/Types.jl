module Types

abstract type AbstractVariable end
abstract type AbstractIntegration end

struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end
struct BoundaryGWForces end
struct BoundaryGWIntegrals end
struct BoundaryGWTendencies end
struct Total end
struct Horizontal end
struct Rho <: AbstractVariable end
struct RhoP <: AbstractVariable end
struct U <: AbstractVariable end
struct V <: AbstractVariable end
struct W <: AbstractVariable end
struct PiP <: AbstractVariable end
struct N2 <: AbstractVariable end
struct DN2DZ <: AbstractVariable end
struct DUDX <: AbstractVariable end
struct DUDY <: AbstractVariable end
struct DUDZ <: AbstractVariable end
struct DVDX <: AbstractVariable end
struct DVDY <: AbstractVariable end
struct DVDZ <: AbstractVariable end
struct Cartesian end
struct TFC end
struct LHS end
struct RHS end
struct EXPL <: AbstractIntegration end
struct IMPL <: AbstractIntegration end
struct X end
struct Y end
struct Z end
struct XZ end
struct YZ end
struct XYZ end

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

include("../MSGWaM/RaySources/activate_orographic_source!.jl")
include("../MSGWaM/RaySources/compute_orographic_mode.jl")
include("../MSGWaM/Interpolation/interpolate_stratification.jl")
include("../MSGWaM/Interpolation/interpolate_mean_flow.jl")
include("../MSGWaM/Interpolation/interpolate.jl")
include("../MSGWaM/Interpolation/get_next_level.jl")
include("../MSGWaM/Interpolation/get_next_half_level.jl")
include("../Update/compute_vertical_wind.jl")
include("../Update/transform.jl")

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

export BoundaryPredictands,
    BoundaryReconstructions,
    BoundaryFluxes,
    BoundaryGWForces,
    BoundaryGWIntegrals,
    BoundaryGWTendencies,
    Total,
    Horizontal,
    Isothermal,
    ConstantCoriolis,
    MCVariant,
    Rho,
    RhoP,
    U,
    V,
    W,
    PiP,
    N2,
    DN2DZ,
    DUDX,
    DUDY,
    DUDZ,
    DVDX,
    DVDY,
    DVDZ,
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
    Cartesian,
    TFC,
    LHS,
    RHS,
    EXPL,
    IMPL,
    X,
    Y,
    Z,
    XZ,
    YZ,
    XYZ,
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
    Rays,
    Increments,
    Integrals,
    SurfaceIndices,
    Forces,
    WKB,
    State

end

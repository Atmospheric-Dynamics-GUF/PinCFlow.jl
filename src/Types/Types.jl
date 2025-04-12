module Types

using MPI

abstract type AbstractBackground end
abstract type AbstractCoriolis end
abstract type AbstractLimiter end
abstract type AbstractVariable end
abstract type AbstractModel end
abstract type AbstractTestCase end
abstract type AbstractBoundaries end
abstract type AbstractSponge end
abstract type AbstractMergeMode end
abstract type AbstractLaunchAlgorithm end
abstract type AbstractIntegration end
abstract type AbstractWKBMode end
abstract type AbstractWKBTestCase <: AbstractTestCase end
abstract type AbstractWKBFilter end

struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end
struct BoundaryGWForces end
struct BoundaryGWIntegrals end
struct BoundaryGWTendencies end

struct Total end
struct Horizontal end

struct Isothermal <: AbstractBackground end

struct ConstantCoriolis <: AbstractCoriolis end

struct MCVariant <: AbstractLimiter end

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

struct PseudoIncompressible <: AbstractModel end

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

struct SteadyState <: AbstractWKBMode end
struct SingleColumn <: AbstractWKBMode end
struct MultiColumn <: AbstractWKBMode end

struct Box <: AbstractWKBFilter end
struct Shapiro <: AbstractWKBFilter end

include("Namelists/DomainNamelist.jl")
include("Namelists/OutputNamelist.jl")
include("Namelists/SettingNamelist.jl")
include("Namelists/DiscretizationNamelist.jl")
include("Namelists/PoissonNamelist.jl")
include("Namelists/AtmosphereNamelist.jl")
include("Namelists/GridNamelist.jl")
include("Namelists/SpongeNamelist.jl")
include("Namelists/WKBNamelist.jl")
include("Namelists/Namelists.jl")

include("Time.jl")
include("Constants.jl")
include("Domain.jl")
include("Grid.jl")
include("Atmosphere.jl")
include("Sponge.jl")

include("Poisson/Tensor.jl")
include("Poisson/Operator.jl")
include("Poisson/Preconditioner.jl")
include("Poisson/BicGStab.jl")
include("Poisson/Correction.jl")
include("Poisson/Poisson.jl")

include("Variables/Predictands.jl")
include("Variables/Tendencies.jl")
include("Variables/Backups.jl")
include("Variables/Auxiliaries.jl")
include("Variables/Reconstructions.jl")
include("Variables/Fluxes.jl")
include("Variables/Variables.jl")

include("WKB/Rays.jl")
include("WKB/Increments.jl")
include("WKB/Integrals.jl")
include("WKB/SurfaceIndices.jl")
include("WKB/Forces.jl")
include("WKB/WKB.jl")

include("State.jl")

include("compute_topography.jl")

# Include functions defined in other modules (needed for initialization). There
# is probably a better way to do this.
include("../Boundaries/set_zonal_boundaries_of_field!.jl")
include("../Boundaries/set_meridional_boundaries_of_field!.jl")
include("../MPIOperations/set_zonal_halos_of_field!.jl")
include("../MPIOperations/set_meridional_halos_of_field!.jl")
include("../MSGWaM/RaySources/activate_orographic_source!.jl")
include("../MSGWaM/RaySources/compute_orographic_mode.jl")
include("../MSGWaM/Interpolation/interpolate_stratification.jl")
include("../MSGWaM/Interpolation/interpolate_mean_flow.jl")

# Export abstract types. Some of these may be moved to other modules.
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

# Export singletons. Some of these may be moved to other modules.
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

# Export composite types.
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
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

struct PseudoIncompressible <: AbstractModel end

struct MountainWave <: AbstractTestCase end
struct WKBMountainWave <: AbstractWKBTestCase end

struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end

struct ExponentialSponge <: AbstractSponge end
struct COSMOSponge <: AbstractSponge end
struct PolynomialSponge <: AbstractSponge end
struct SinusoidalSponge <: AbstractSponge end

struct ConstantWaveEnergy <: AbstractMergeMode end
struct ConstantWaveAction <: AbstractMergeMode end

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

struct SteadyState <: AbstractWKBMode end
struct SingleColumn <: AbstractWKBMode end
struct MultiColumn <: AbstractWKBMode end

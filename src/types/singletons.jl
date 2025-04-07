struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end

struct Total end
struct Horizontal end

struct Isothermal <: AbstractBackground end

struct ConstantCoriolis <: AbstractCoriolis end

struct MCVariant <: AbstractLimiter end

struct Rho <: AbstractVariable end
struct RhoP <: AbstractVariable end
struct U <: AbstractVariable end
struct US <: AbstractVariable end
struct V <: AbstractVariable end
struct VS <: AbstractVariable end
struct W <: AbstractVariable end
struct WS <: AbstractVariable end
struct WTFC <: AbstractVariable end
struct WSTFC <: AbstractVariable end
struct ThetaP <: AbstractVariable end
struct PiP <: AbstractVariable end

struct PseudoIncompressible <: AbstractModel end

struct MountainWave <: AbstractTestCase end
struct RayTracer <: AbstractTestCase end

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

struct xDir <: AbstractDir end
struct yDir <: AbstractDir end
struct zDir <: AbstractDir end
abstract type AbstractVariable end
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

abstract type AbstractCoordinate end
struct Cartesian <: AbstractCoordinate end
struct TFC <: AbstractCoordinate end
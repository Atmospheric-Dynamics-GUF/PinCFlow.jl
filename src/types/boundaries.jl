abstract type AbstractBoundaries end
struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end

abstract type AbstractBoundaryVariables end
struct BoundaryPredictands <: AbstractBoundaryVariables end
struct BoundaryReconstructions <: AbstractBoundaryVariables end
struct BoundaryFluxes <: AbstractBoundaryVariables end
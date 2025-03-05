using Trixi: AbstractEquations

abstract type AbstractMatrixSolver end
abstract type AbstractSpatialSolver end
abstract type AbstractBoundaryCondition end

struct FiniteVolumeSolver <: AbstractSpatialSolver end

struct PeriodicBC <: AbstractBoundaryCondition end
struct SolidWallBC <: AbstractBoundaryCondition end

struct BoundaryConditions{LeftBC, RightBC, BottomBC, TopBC}
    left::LeftBC
    right::RightBC
    bottom::BottomBC
    top::TopBC
    function BoundaryConditions(left, right, bottom, top)
        if left isa PeriodicBC ||
           right isa PeriodicBC ||
           bottom isa PeriodicBC ||
           top isa PeriodicBC
            @assert left isa PeriodicBC &&
                    right isa PeriodicBC &&
                    bottom isa PeriodicBC &&
                    top isa PeriodicBC
        end
        new{typeof(left), typeof(right), typeof(bottom), typeof(top)}(left, right, bottom,
                                                                      top)
    end
end

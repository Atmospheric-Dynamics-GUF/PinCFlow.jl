struct SemiDiscretization{
    Grid,
    Equations<:AbstractEquations,
    SurfaceFlux,
    IC,
    BC,
    Solver,
    MatrixSolver<:AbstractMatrixSolver,
    Cache,
}
    grid::Grid
    equations::Equations
    surface_flux::SurfaceFlux
    initial_condition::IC
    boundary_conditions::BC
    solver::Solver
    matrix_solver::MatrixSolver
    cache::Cache
end

function SemiDiscretization(
    grid,
    equations,
    surface_flux,
    initial_condition;
    solver = FiniteVolumeSolver(),
    boundary_conditions = BoundaryConditions(
        PeriodicBC(),
        PeriodicBC(),
        PeriodicBC(),
        PeriodicBC(),
    ),
    matrix_solver = CGSolver(),
    cache = (;),
)
    cache = (; cache..., create_cache(equations, grid, initial_condition)...)

    SemiDiscretization(
        grid,
        equations,
        surface_flux,
        initial_condition,
        boundary_conditions,
        solver,
        matrix_solver,
        cache,
    )
end

function create_cache(equations, grid, initial_condition)
    # TODO - Add stuff here
    cache = (;)

    return cache
end

"""
```julia
apply_corrector!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
```

Apply the pressure correction step to enforce mass conservation.

This is the main interface for the pressure correction procedure. It computes the
right-hand side, solves the Poisson equation, and applies velocity/density corrections
to satisfy the divergence-free condition.

# Arguments

  - `state`: Complete simulation state
  - `dt`: Time step size
  - `facray`: Rayleigh damping factor for sponge boundaries
  - `facprs`: Pressure correction scaling factor

# Returns

  - `(errflagbicg, niterbicg)`: Error status and iteration count from linear solver

# Process

 1. Initialize and compute RHS of Poisson equation
 2. Solve linear system using BiCGStab
 3. Set boundary conditions on pressure correction
 4. Apply momentum and buoyancy corrections

# See also

  - [`PinCFlow.PoissonSolver.compute_rhs!`](@ref)
  - [`PinCFlow.PoissonSolver.solve_poisson!`](@ref)
  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
  - [`PinCFlow.PoissonSolver.correct!`](@ref)
"""
function apply_corrector!(
    state::State,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; namelists, domain) = state
    (; model, zboundaries) = namelists.setting
    (; rhs) = state.poisson
    (; dpip) = state.variables.tendencies

    # Initialize RHS.
    rhs .= 0.0

    # Calculate RHS and tolreance reference.
    tolref = compute_rhs!(state, rhs, model)

    # Solve Poisson equation.
    (errflagbicg, niterbicg) =
        solve_poisson!(state, rhs, tolref, dt, facray, facprs)

    # Return if an error occurred.
    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    # Set boundaries of pressure correction.
    set_zonal_boundaries_of_field!(dpip, namelists, domain)
    set_meridional_boundaries_of_field!(dpip, namelists, domain)
    set_vertical_boundaries_of_field!(dpip, namelists, domain, zboundaries, +)

    # Correct momentum and buoyancy.
    correct!(state, dt, facray, facprs)

    return (errflagbicg, niterbicg)
end

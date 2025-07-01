"""
```julia
State{
    A <: Namelists,
    B <: Time,
    C <: Constants,
    D <: Domain,
    E <: Grid,
    F <: Atmosphere,
    G <: Sponge,
    H <: Poisson,
    I <: Variables,
    J <: WKB,
}
```

Main simulation state container.

Holds all components required for the simulation including physical
parameters, computational grids, field variables, and specialized solvers.

# Type Parameters

  - `A<:Namelists`: Configuration and parameter specifications
  - `B<:Time`: Time integration scheme parameters
  - `C<:Constants`: Physical constants and reference scales
  - `D<:Domain`: MPI domain decomposition and grid indices
  - `E<:Grid`: Computational grid and coordinate system
  - `F<:Atmosphere`: Background atmospheric state
  - `G<:Sponge`: Boundary sponge layer configuration
  - `H<:Poisson`: Pressure correction solver components
  - `I<:Variables`: All prognostic and diagnostic field variables
  - `J<:WKB`: Gravity wave parameterization (WKB/ray tracing)

# Fields

  - `namelists::A`: Simulation configuration and physical parameters
  - `time::B`: Runge-Kutta time integration coefficients
  - `constants::C`: Physical constants and non-dimensionalization scales
  - `domain::D`: Parallel domain decomposition and MPI communication
  - `grid::E`: Computational mesh, coordinates, and metric tensors
  - `atmosphere::F`: Background state (pressure, density, stratification)
  - `sponge::G`: Absorbing boundary layer damping coefficients
  - `poisson::H`: Pressure Poisson solver (BiCGStab, preconditioner, operators)
  - `variables::I`: All field variables (velocities, pressure, density, etc.)
  - `wkb::J`: Gravity wave ray tracing and parameterization

# Constructor

```julia
State(namelists::Namelists)
```

Initialize complete simulation state from configuration parameters.

Creates and links all simulation components in proper dependency order:

 1. Physical constants and reference scales
 2. Time integration parameters
 3. MPI domain decomposition
 4. Computational grid and coordinates
 5. Background atmospheric state
 6. Boundary conditions and damping
 7. Linear solvers and operators
 8. Field variable storage
 9. Parameterization schemes

# Usage Context

## Initialization

```julia
# Load configuration
namelists = Namelists("config.toml")

# Initialize complete simulation state
state = State(namelists)

# Access components
grid = state.grid
variables = state.variables
constants = state.constants
```

## Field Access

```julia
# Prognostic variables
(; u, v, w) = state.variables.predictands
(; rho, p) = state.variables.predictands

# Background state
p0 = state.atmosphere.pstrattfc
rho0 = state.atmosphere.rhostrattfc

# Grid information
(; dx, dy, dz) = state.grid
```

# Component Relationships

## Data Flow

  - `namelists` → configures all other components
  - `constants` → sets reference scales for non-dimensionalization
  - `domain` → defines parallel decomposition for all fields
  - `grid` → provides coordinate system for all spatial operations
  - `atmosphere` → supplies background state for perturbation dynamics
  - `variables` → stores evolving solution fields
  - `poisson` → solves pressure correction using grid/domain info
  - `sponge` → applies boundary damping to variables
  - `wkb` → parameterizes subgrid gravity waves

## Memory Layout

  - All field arrays sized according to `domain` specifications
  - Grid coordinates computed from `namelists.domain` parameters
  - Background profiles initialized from `namelists.atmosphere`
  - Solver workspace allocated based on `namelists.poisson` settings

## Parallel Computing

  - `domain` handles MPI communication patterns
  - All field operations respect domain decomposition
  - Boundary exchanges coordinated through `domain.comm`
  - Solver convergence computed via global reductions

# Notes

  - Single point of access for entire simulation state
  - Ensures consistent initialization order and dependencies
  - All components initialized from single configuration source

# See also

  - [`Namelists`](@ref)
  - [`Variables`](@ref)
  - [`Grid`](@ref)
  - [`Poisson`](@ref)
"""
struct State{
    A <: Namelists,
    B <: Time,
    C <: Constants,
    D <: Domain,
    E <: Grid,
    F <: Atmosphere,
    G <: Sponge,
    H <: Poisson,
    I <: Variables,
    J <: WKB,
}
    namelists::A
    time::B
    constants::C
    domain::D
    grid::E
    atmosphere::F
    sponge::G
    poisson::H
    variables::I
    wkb::J
end

"""
```julia
State(namelists::Namelists)
```

Initialize complete simulation state from configuration.

Constructs all simulation components in dependency order, ensuring proper
initialization of physical constants, computational grid, background state,
field storage, and numerical solvers.

# Arguments

  - `namelists::Namelists`: Complete simulation configuration

# Returns

  - `State`: Fully initialized simulation state ready for time integration

# Initialization Order

 1. **Constants**: Physical parameters and reference scales
 2. **Time**: Runge-Kutta integration coefficients
 3. **Domain**: MPI decomposition and communication setup
 4. **Grid**: Computational mesh and coordinate transformations
 5. **Atmosphere**: Background state profiles and stratification
 6. **Sponge**: Boundary layer damping configuration
 7. **Poisson**: Pressure solver setup (operators, preconditioner, workspace)
 8. **Variables**: Field variable storage (predictands, tendencies, auxiliaries)
 9. **WKB**: Gravity wave parameterization initialization

# Example

```julia
# Initialize from configuration file
namelists = Namelists("simulation.toml")
state = State(namelists)

# Ready for time integration
integrate!(state, total_time)
```
"""
function State(namelists::Namelists)

    # Initialize everything.
    constants = Constants(namelists)
    time = Time()
    domain = Domain(namelists)
    grid = Grid(namelists, constants, domain)
    atmosphere = Atmosphere(namelists, constants, domain, grid)
    sponge = Sponge(namelists, domain, grid)
    poisson = Poisson(namelists, domain)
    variables = Variables(namelists, constants, domain, atmosphere)
    wkb = WKB(namelists, constants, domain, grid)

    # Return a State instance.
    return State(
        namelists,
        time,
        constants,
        domain,
        grid,
        atmosphere,
        sponge,
        poisson,
        variables,
        wkb,
    )
end

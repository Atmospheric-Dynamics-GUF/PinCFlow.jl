"""
```julia
Variables{
    A <: Predictands,
    B <: Tendencies,
    C <: Backups,
    D <: Auxiliaries,
    E <: Reconstructions,
    F <: Fluxes,
}
```

Complete field variable storage for the simulation.

Central container for all field variables and computational arrays used during
time integration, including prognostic fields, tendencies, auxiliary variables,
and flux calculations.

# Type Parameters

  - `A<:Predictands`: Prognostic variable storage
  - `B<:Tendencies`: Time tendency storage
  - `C<:Backups`: Backup copies for time stepping
  - `D<:Auxiliaries`: Diagnostic and auxiliary fields
  - `E<:Reconstructions`: MUSCL reconstructed interface values
  - `F<:Fluxes`: Flux arrays for finite volume calculations

# Fields

  - `predictands::A`: Prognostic variables (u, v, w, ρ, p, π')
  - `tendencies::B`: Time derivatives for Runge-Kutta integration
  - `backups::C`: Previous time step storage for multi-stage scheme
  - `auxiliaries::D`: Diagnostic fields and intermediate calculations
  - `reconstructions::E`: Left/right states at cell interfaces (MUSCL)
  - `fluxes::F`: Flux vectors for advective and pressure terms

# Constructor

    Variables(namelists, constants, domain, atmosphere)

Create complete variable storage from simulation configuration.

# Usage Context

```julia
# Initialize all field storage
vars = Variables(namelists, constants, domain, atmosphere)

# Access prognostic fields
u, v, w = vars.predictands.u, vars.predictands.v, vars.predictands.w
ρ, p = vars.predictands.rho, vars.predictands.p

# Access tendencies during RK stages
du_dt = vars.tendencies.du
dp_dt = vars.tendencies.dpip

# Use reconstructions in flux calculations
apply_3d_muscl!(vars.predictands, vars.reconstructions, state)
compute_fluxes!(vars.fluxes, vars.reconstructions, state)
```

# Component Relationships

  - `predictands` stores evolving solution fields
  - `tendencies` accumulates time derivative contributions
  - `backups` preserves previous values for multi-stage time stepping
  - `auxiliaries` computes diagnostic quantities from predictands
  - `reconstructions` provides high-order interface values for fluxes
  - `fluxes` stores numerical flux vectors for finite volume updates

# Memory Organization

  - All arrays sized according to domain decomposition
  - Includes halo regions for MPI boundary exchanges
  - Model-dependent allocation (Boussinesq vs Compressible)
"""
struct Variables{
    A <: Predictands,
    B <: Tendencies,
    C <: Backups,
    D <: Auxiliaries,
    E <: Reconstructions,
    F <: Fluxes,
}
    predictands::A
    tendencies::B
    backups::C
    auxiliaries::D
    reconstructions::E
    fluxes::F
end

"""
```julia
Variables(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)
```

Initialize complete field variable storage for atmospheric simulation.

Creates and organizes all variable containers required for the simulation,
including prognostic fields, time derivatives, auxiliary arrays, and
computational workspace for finite volume calculations.

# Arguments

  - `namelists::Namelists`: Simulation configuration
  - `constants::Constants`: Physical constants and reference scales
  - `domain::Domain`: Computational domain and MPI decomposition
  - `atmosphere::Atmosphere`: Background atmospheric state

# Returns

  - `Variables`: Complete field variable storage
"""
function Variables(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)

    # Initialize all fields.
    predictands = Predictands(namelists, constants, domain, atmosphere)
    tendencies = Tendencies(namelists, domain)
    backups = Backups(domain)
    auxiliaries = Auxiliaries(domain)
    reconstructions = Reconstructions(domain)
    fluxes = Fluxes(namelists, domain)

    # Return a Variables instance.
    return Variables(
        predictands,
        tendencies,
        backups,
        auxiliaries,
        reconstructions,
        fluxes,
    )
end

"""
    set_boundary_rays!(state::State)

Entry point for setting ray boundary conditions based on test case type.

Dispatches to the appropriate boundary condition implementation depending on
whether the simulation uses WKB ray tracing.

# Arguments

  - `state::State`: Complete simulation state containing ray data
"""
function set_boundary_rays!(state::State)
    (; testcase) = state.namelists.setting
    set_boundary_rays!(state, testcase)
    return
end

"""
    set_boundary_rays!(state::State, testcase::AbstractTestCase)

No-op for non-WKB test cases.

Standard test cases don't use ray tracing, so no ray boundary conditions are needed.

# Arguments

  - `state::State`: Simulation state (unused)
  - `testcase::AbstractTestCase`: Non-WKB test case
"""
function set_boundary_rays!(state::State, testcase::AbstractTestCase)
    return
end

"""
    set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)

Set ray boundary conditions for WKB test cases based on WKB mode.

Dispatches to the specific WKB mode implementation for ray boundary handling.

# Arguments

  - `state::State`: Simulation state containing WKB configuration
  - `testcase::AbstractWKBTestCase`: WKB test case specification
"""
function set_boundary_rays!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb
    set_boundary_rays!(state, wkb_mode)
    return
end

"""
    set_boundary_rays!(state::State, wkb_mode::SteadyState)

Set ray boundary conditions for steady-state WKB mode.

Applies horizontal boundary conditions only, as steady-state mode typically
handles vertical boundaries differently (no reflection/absorption needed).

# Arguments

  - `state::State`: Simulation state
  - `wkb_mode::SteadyState`: Steady-state WKB mode

# Boundary Conditions Applied

  - **Zonal**: Periodic or wall boundaries in x-direction (if sizex > 1)
  - **Meridional**: Periodic or wall boundaries in y-direction (if sizey > 1)
  - **Vertical**: No special handling (rays transported vertically)

# Steady-State Assumptions

  - Background atmosphere is time-independent
  - Vertical boundary effects handled by propagation scheme
  - Only horizontal domain coupling needed for MPI decomposition
"""
function set_boundary_rays!(state::State, wkb_mode::SteadyState)
    (; sizex, sizey) = state.namelists.domain

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end

    return
end

"""
    set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)

Set ray boundary conditions for general WKB modes.

Applies comprehensive boundary conditions in all active spatial dimensions,
including vertical boundaries with reflection, absorption, or transmission.

# Arguments

  - `state::State`: Simulation state
  - `wkb_mode::AbstractWKBMode`: General WKB mode (MultiColumn, SingleColumn, etc.)

# Boundary Conditions Applied

  - **Zonal**: Domain boundaries in x-direction (if sizex > 1)
  - **Meridional**: Domain boundaries in y-direction (if sizey > 1)
  - **Vertical**: Physical boundaries (surface reflection, top absorption)

# Vertical Boundary Types

Determined by `zboundaries` setting:

  - **SolidWallBoundaries**: Reflection at surface, absorption at top
  - **PeriodicBoundaries**: Vertical periodicity (rare, for idealized cases)
  - **OpenBoundaries**: Transmission through boundaries

# Physical Implementation

## Surface Reflection

  - Rays hitting topography are reflected upward
  - Vertical wavenumber sign flipped: `m → -m`
  - Conserves wave energy while changing propagation direction

## Top Absorption

  - Rays reaching domain top are removed or damped
  - Prevents artificial wave accumulation
  - Simulates radiative damping in upper atmosphere

# MPI Coordination

Ensures proper ray exchange between:

  - Horizontally adjacent processes
  - Vertically stacked processes
  - Corner/edge processes in 2D/3D decompositions

# Order of Operations

 1. Apply zonal boundaries (x-direction)
 2. Apply meridional boundaries (y-direction)
 3. Apply vertical boundaries (z-direction)

This order ensures all domain coupling is properly handled.
"""
function set_boundary_rays!(state::State, wkb_mode::AbstractWKBMode)
    (; sizex, sizey) = state.namelists.domain
    (; zboundaries) = state.namelists.setting

    if sizex > 1
        set_zonal_boundary_rays!(state)
    end
    if sizey > 1
        set_meridional_boundary_rays!(state)
    end
    set_vertical_boundary_rays!(state, zboundaries)

    return
end

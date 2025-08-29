"""
```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
)
```

Enforce vertical boundary conditions for all predictand fields.

The symmetry conditions are as follows:

  - Density-fluctuations fields (`rho`, `rhop`): point reflection (`-`)

  - Vertical velocity (`w`): point reflection (`-`) on the staggered grid

  - Horizontal velocities (`u`, `v`): line reflection (`+`)

  - Exner-pressure fluctuations (`pip`): line reflection (`+`)

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
```

Enforce vertical boundary conditions for all reconstruction fields.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
)
```

Set the vertical fluxes at the vertical boundaries to zero.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
)
```

Enforce vertical boundary conditions for gravity-wave-integral fields by dispatching to a WKB-mode-specific method.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
)
```

Enforce vertical boundary conditions for gravity-wave-integral fields needed in `SingleColumn` and `SteadyState` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::MultiColumn,
)
```

Enforce vertical boundary conditions for gravity-wave-integral fields needed in `MultiColumn` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
)
```

Enforce vertical boundary conditions for gravity-wave-tendency fields by dispatching to a WKB-mode-specific method.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
)
```

Enforce vertical boundary conditions for gravity-wave-tendency fields needed in `SingleColumn` and `SteadyState` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::MultiColumn,
)
```

Enforce vertical boundary conditions for gravity-wave-tendency fields needed in `MultiColumn` configurations, using line reflection.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_compressible_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_tracer_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_ice_vertical_boundaries!`](@ref)

  - [`PinCFlow.Boundaries.set_turbulence_vertical_boundaries!`](@ref)
"""
function set_vertical_boundaries! end

function set_vertical_boundaries!(state::State, variables::BoundaryPredictands)
    (; namelists, domain) = state
    (; model) = namelists.setting
    (; rho, rhop, u, v, w, pip) = state.variables.predictands
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

    set_vertical_boundaries_of_field!(rho, namelists, domain, -)
    set_vertical_boundaries_of_field!(rhop, namelists, domain, -)

    set_vertical_boundaries_of_field!(w, namelists, domain, -; staggered = true)

    set_vertical_boundaries_of_field!(u, namelists, domain, +)
    set_vertical_boundaries_of_field!(v, namelists, domain, +)

    set_vertical_boundaries_of_field!(pip, namelists, domain, +)

    set_compressible_vertical_boundaries!(state, variables, model)

    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

    for field in fieldnames(Reconstructions)
        set_vertical_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

function set_vertical_boundaries!(state::State, variables::BoundaryFluxes)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables
    (; model) = state.namelists.setting
    (; tracersetup) = state.namelists.tracer
    (; icesetup) = state.namelists.ice
    (; turbulencesetup) = state.namelists.turbulence

    # Set all vertical boundary fluxes to zero.

    if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    if ko + nzz == sizezz
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    set_compressible_vertical_boundaries!(state, variables, model)
    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

function set_vertical_boundaries!(state::State, variables::BoundaryWKBIntegrals)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, wkb_mode)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_vertical_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain,
            +;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_vertical_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain,
            +;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, wkb_mode)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_vertical_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
            +,
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_vertical_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
            +,
        )
    end

    return
end

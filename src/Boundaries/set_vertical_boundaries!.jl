"""
```julia
set_vertical_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        BoundaryFluxes,
    },
)
```

Enforce vertical boundary conditions for predictands, reconstructions or fluxes, by dispatching to the appropriate method.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Boussinesq,
)
```

Enforce vertical boundary conditions for predictands in Boussinesq mode.

The symmetry conditions are as follows:

  - Density fluctuations (`rhop`): point reflection (`-`)

  - Vertical velocity (`w`): point reflection (`-`) on the staggered grid

  - Horizontal velocities (`u`, `v`): line reflection (`+`)

  - Exner-pressure fluctuations (`pip`): line reflection (`+`)

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::PseudoIncompressible,
)
```

Enforce vertical boundary conditions for predictands in pseudo-incompressible mode.

In contrast to Boussinesq mode, this includes the density (`rho`). Since it is offset by the background, it is point-reflected, like the density fluctuations.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
```

Enforce vertical boundary conditions for predictands in compressible mode.

In contrast to pseudo-incompressible modes, this includes the mass-weighted potential temperature (`p`), which is line-reflected.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Boussinesq,
)
```

Enforce vertical boundary conditions for reconstructions in Boussinesq mode.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{PseudoIncompressible, Compressible},
)
```

Enforce vertical boundary conditions for reconstructions in non-Boussinesq modes.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Boussinesq,
)
```

Enforce vertical boundary conditions for vertical fluxes in Boussinesq mode.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::PseudoIncompressible,
)
```

Enforce vertical boundary conditions for vertical fluxes in pseudo-incompressible mode.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
```

Enforce vertical boundary conditions for vertical fluxes in compressible mode.

```julia
set_vertical_boundaries!(state::State, variables::AbstractBoundaryWKBVariables)
```

Enforce vertical boundary conditions for WKB variables by dispatching to the appropriate method.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn},
)
```

Enforce vertical boundary conditions for WKB integrals needed in `SingleColumn` and `SteadyState` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::MultiColumn,
)
```

Enforce vertical boundary conditions for WKB integrals needed in `MultiColumn` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::Union{SteadyState, SingleColumn},
)
```

Enforce vertical boundary conditions for WKB tendencies needed in `SingleColumn` and `SteadyState` configurations, using line reflection.

```julia
set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBTendencies,
    wkb_mode::MultiColumn,
)
```

Enforce vertical boundary conditions for WKB tendencies needed in `MultiColumn` configurations, using line reflection.

# Arguments

  - `state`: Model state.

  - `variables`: Boundary-variable category.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
"""
function set_vertical_boundaries! end

function set_vertical_boundaries!(
    state::State,
    variables::Union{
        BoundaryPredictands,
        BoundaryReconstructions,
        BoundaryFluxes,
    },
)
    (; model) = state.namelists.atmosphere
    set_vertical_boundaries!(state, variables, model)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Boussinesq,
)
    (; namelists, domain) = state
    (; rhop, u, v, w, pip) = state.variables.predictands

    set_vertical_boundaries_of_field!(rhop, namelists, domain, -)
    set_vertical_boundaries_of_field!(u, namelists, domain, +)
    set_vertical_boundaries_of_field!(v, namelists, domain, +)
    set_vertical_boundaries_of_field!(w, namelists, domain, -; staggered = true)
    set_vertical_boundaries_of_field!(pip, namelists, domain, +)

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::PseudoIncompressible,
)
    (; namelists, domain) = state
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    set_vertical_boundaries_of_field!(rho, namelists, domain, -)
    set_vertical_boundaries_of_field!(rhop, namelists, domain, -)
    set_vertical_boundaries_of_field!(u, namelists, domain, +)
    set_vertical_boundaries_of_field!(v, namelists, domain, +)
    set_vertical_boundaries_of_field!(w, namelists, domain, -; staggered = true)
    set_vertical_boundaries_of_field!(pip, namelists, domain, +)

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands

    set_vertical_boundaries_of_field!(rho, namelists, domain, -)
    set_vertical_boundaries_of_field!(rhop, namelists, domain, -)
    set_vertical_boundaries_of_field!(u, namelists, domain, +)
    set_vertical_boundaries_of_field!(v, namelists, domain, +)
    set_vertical_boundaries_of_field!(w, namelists, domain, -; staggered = true)
    set_vertical_boundaries_of_field!(pip, namelists, domain, +)
    set_vertical_boundaries_of_field!(p, namelists, domain, +)

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Boussinesq,
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in (:rhoptilde, :utilde, :vtilde, :wtilde)
        set_vertical_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    model::Union{PseudoIncompressible, Compressible},
)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_vertical_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Boussinesq,
)
    (; z_size) = state.namelists.domain
    (; nz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables

    @ivy if ko == 0
        for field in (:phirhop, :phiu, :phiv, :phitheta)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    @ivy if ko + nz == z_size
        for field in (:phirhop, :phiu, :phiv, :phiw, :phitheta)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::PseudoIncompressible,
)
    (; z_size) = state.namelists.domain
    (; nz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables

    @ivy if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv, :phitheta)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    @ivy if ko + nz == z_size
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw, :phitheta)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
    (; z_size) = state.namelists.domain
    (; nz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables

    @ivy if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv, :phitheta, :phip)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    @ivy if ko + nz == z_size
        for field in fieldnames(Fluxes)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryWKBVariables,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, wkb_mode)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryWKBIntegrals,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e, :sterm, :bterm)
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

    for field in (:uu, :uv, :uw, :vv, :vw, :utheta, :vtheta, :e, :sterm, :bterm)
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
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dtkedt)
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

    for field in (:dudt, :dvdt, :dthetadt, :dtkedt)
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
    variables::BoundaryDiffusionCoefficients,
)
    (; namelists, domain) = state
    (; turbulencediffusioncoefficients) = state.turbulence

    for field in (:kh, :km, :kek)
        set_vertical_boundaries_of_field!(
            getfield(turbulencediffusioncoefficients, field),
            namelists,
            domain,
            -,
        )
    end

    return
end

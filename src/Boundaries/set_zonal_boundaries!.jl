"""
```julia
set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
```

Enforce zonal boundary conditions for all predictand fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
  - [`PinCFlow.Boundaries.set_compressible_zonal_boundaries!`](@ref)
"""
function set_zonal_boundaries!(state::State, variables::BoundaryPredictands)
    (; namelists, domain) = state
    (; predictands) = state.variables
    (; model) = namelists.setting

    for field in (:rho, :rhop, :u, :v, :w, :pip)
        set_zonal_boundaries_of_field!(
            getfield(predictands, field),
            namelists,
            domain,
        )
    end

    set_compressible_zonal_boundaries!(state, model)

    return
end

"""
```julia
set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
```

Enforce zonal boundary conditions for all reconstruction fields.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries!(state::State, variables::BoundaryReconstructions)
    (; namelists, domain) = state
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_zonal_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

"""
```julia
set_zonal_boundaries!(state::State, variables::BoundaryGWIntegrals)
```

Enforce zonal boundary conditions for gravity-wave-integral fields, dispatching based on WKB mode.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryGWIntegrals)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

"""
```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
)
```

Enforce zonal boundary conditions for gravity-wave-integral fields needed in `SingleColumn` and `SteadyState` configurations.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
```

Enforce zonal boundary conditions for gravity-wave-integral fields needed in `MultiColumn` configurations.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_zonal_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
```julia
set_zonal_boundaries!(state::State, variables::BoundaryGWTendencies)
```

Enforce zonal boundary conditions for gravity-wave-tendency fields, dispatching based on WKB mode.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
"""
function set_zonal_boundaries!(state::State, variables::BoundaryGWTendencies)
    (; wkb_mode) = state.namelists.wkb
    set_zonal_boundaries!(state, variables, wkb_mode)
    return
end

"""
```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
)
```

Enforce zonal boundary conditions for gravity-wave-tendency fields needed in `SingleColumn` and `SteadyState` configurations.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

"""
```julia
set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::MultiColumn,
)
```

Enforce zonal boundary conditions for gravity-wave-tendency fields needed in `MultiColumn` configurations.

# Arguments

  - `state`: Model state.
  - `variables`: Boundary-variable category.
  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
"""
function set_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_zonal_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
        )
    end

    return
end

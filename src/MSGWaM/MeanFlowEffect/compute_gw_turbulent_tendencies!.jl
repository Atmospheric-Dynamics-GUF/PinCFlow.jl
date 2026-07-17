"""
```julia
compute_gw_turbulent_tendencies!(state::State, i::Integer, j::Integer, k::Integer)
```

Calculates the tendency that is to be added to the turbulence equation, by dispatching to the appropriate method.

```julia
compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
```

Compute and return the gravity-wave turbulent tendencies at ``\\left(i, j, k\\right)``.

Calculates the turbulent tendency from the gravity-wave induced shear as follows,

```math
\\begin{align*}
    \\left(\\frac{\\partial e_\\mathrm{k}}{\\partial t}\\right)_\\mathrm{w} = K_\\mathrm{M} \\ \\mathcal{S}_{gw} \\;,
\\end{align*}
```

where ``K_\\mathrm{H}`` represents the eddy diffusion coefficient for momentum. (see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_shear!`](@ref) for documentation on the gravity-wave induced shear).

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `turbulence_scheme`: General turbulence parameterization configuration.

!!! danger "Experimental"
    The gravity-wave induced shear is an experimental feature that hasn't been validated yet.
"""
function compute_gw_turbulent_tendencies! end

function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulent_tendencies!(
        state,
        i,
        j,
        k,
        Val(turbulence_scheme),
    )
    return
end

function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

@ivy function compute_gw_turbulent_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; gw_shear) = state.turbulence.turbulencewkbtendencies

    gw_shear[i, j, k] = turbulence_diffusion_coefficient(state, i, j, k, KM()) * 
        gw_shear[i, j, k]

    return
end

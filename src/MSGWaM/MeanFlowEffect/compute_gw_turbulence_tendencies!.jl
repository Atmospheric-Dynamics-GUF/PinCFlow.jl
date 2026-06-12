"""
```julia
compute_gw_turbulence_tendencies!(state::State, i::Integer, j::Integer, k::Integer)
```

Compute the turbulence production due to gravity wave shear by dispatching to the appropriate method.

```julia
compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
```

Return for configurations without turbulence parameterization.

```julia
compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
```

Compute and return the turbulence production due to gravity wave shear at ``\\left(i, j, k\\right)``.

Calculates the tendency that is to be added to the TKE equation, given by

```math
\\begin{align*}
    \\left(\\frac{\\partial \\rho_\\mathrm{b} e_\\mathrm{k,b}}{\\partial t}\\right)_\\mathrm{w} & = \\rho_\\mathrm{b}K_{e_\\mathrm{k}}S_\\mathrm{w} ,
\\end{align*}
```

where ``e_\\mathrm{k,b}`` is the resolved TKE and ``\\rho_\\mathrm{b}`` is the resolved density (including the reference part ``\\bar{\\rho}``). For a documentation of the gravity wave shear ``S_\\mathrm{w}``, see [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_turbulence_integrals!`](@ref).

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `turbulence_scheme`: General turbulence-transport configuration.
"""
function compute_gw_turbulence_tendencies! end

function compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulence_tendencies!(
        state,
        i,
        j,
        k,
        Val(turbulence_scheme),
    )
    return
end

function compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

@ivy function compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; gwshear, gwbuoy) = state.turbulence.turbulencewkbintegrals
    (; dtkedt) = state.turbulence.turbulencewkbtendencies
    (; wave_impact) = state.namelists.turbulence

    if !wave_impact
        return
    end

    dtkedt[i, j, k] =
        turbulence_diffusion_coefficient(state, i, j, k, KM()) *
        gwshear[i, j, k]
    # + turbulence_diffusion_coefficient(state, i, j, k, KH()) * gwbuoy[i, j, k]

    return
end

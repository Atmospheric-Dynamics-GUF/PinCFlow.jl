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

function compute_gw_turbulence_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; dtkedt) = state.turbulence.turbulencewkbtendencies
    (; shear) = state.turbulence.turbulencewkbintegrals
    (; km) = state.turbulence.turbulencediffusioncoefficients
    (; rhobar) = state.atmosphere

    dtkedt[i, j, k] = shear[i, j, k]  # km[i, j, k] * shear[i, j, k] #/ rhobar[i, j, k]

    return
end
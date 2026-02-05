function compute_turbulence_forcing! end

function compute_turbulence_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    compute_turbulence_forcing!(state, i, j, k, turbulence_scheme)

    return
end

function compute_turbulence_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::NoTurbulence,
)
    return
end

function compute_turbulence_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::TKEScheme,
)
    (; km, kh) = state.turbulence.turbulencediffusioncoefficients
    (; dtkedt) = state.wkb.tendencies
    (; sterm, bterm) = state.wkb.integrals
    
    @ivy dtkedt[i, j, k] =
        km[i, j, k] * sterm[i, j, k] + kh[i, j, k] * bterm[i, j, k]

    return
end
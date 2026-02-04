function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
    (; turbulence_scheme) = state.namelists.turbulence

    return compute_turbulent_damping(state, r, i, j, k, turbulence_scheme)
end

function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::NoTurbulence,
)::AbstractFloat
    return 0.0
end

function compute_turbulent_damping(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::TKEScheme,
)::AbstractFloat

    (; rays) = state.wkb 

    
    return 0.0
end
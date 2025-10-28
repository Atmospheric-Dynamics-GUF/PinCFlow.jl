function compute_diffusion_coefficients end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
)::AbstractFloat
    (; kinematic_diffusivity) = state.namelists.atmosphere
    (; uref, lref) = state.constants
    (; rhobar) = state.atmosphere
    (; k0) = state.domain
    (; turbulence_scheme) = state.namelists.turbulence

    mu_mom_diff = kinematic_diffusivity / uref / lref
    coef_d = mu_mom_diff * rhobar[i, j, k0]
    coef_d += compute_diffusion_coefficients(
        state,
        i,
        j,
        k,
        variable,
        turbulence_scheme,
    )

    return coef_d
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    k::Integer,
    variable::Union{U, V, W, Theta, Chi},
    turbulence_scheme::NoTurbulence,
)::AbstractFloat
    return 0.0
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
    turbulence_scheme::TKEScheme,
)::AbstractFloat
    (; km) = state.turbulence.turbulencediffusioncoefficients
    (; momentum_coupling) = state.namelists.turbulence

    if momentum_coupling
        return km[i, j, k]
    else
        return 0.0
    end
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Theta,
    turbulence_scheme::TKEScheme,
)::AbstractFloat
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; entropy_coupling) = state.namelists.turbulence

    if entropy_coupling
        return kh[i, j, k]
    else
        return 0.0
    end
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Chi,
)::AbstractFloat
    (; turbulence_scheme) = state.namelists.turbulence

    return compute_diffusion_coefficients(
        state,
        i,
        j,
        k,
        variable,
        turbulence_scheme,
    )
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Chi,
    turbulence_scheme::TKEScheme,
)::AbstractFloat
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    return kh[i, j, k]
end

function compute_diffusion_coefficients(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Theta,
)::AbstractFloat
    (; thermal_conductivity) = state.namelists.atmosphere
    (; uref, lref) = state.constants
    (; rhobar) = state.atmosphere
    (; k0) = state.domain
    (; turbulence_scheme) = state.namelists.turbulence

    mu_conduct = thermal_conductivity / uref / lref
    coef_t = mu_conduct * rhobar[i, j, k0] / rhobar[i, j, k]
    coef_t += compute_diffusion_coefficients(
        state,
        i,
        j,
        k,
        variable,
        turbulence_scheme,
    )

    return coef_t
end
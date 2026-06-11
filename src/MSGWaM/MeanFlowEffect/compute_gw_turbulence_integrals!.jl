"""
```julia
compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the gravity wave shear by dispatching to the appropriate method.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:NoTurbulence},
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Return for configurations without turbulence scheme.

```julia
compute_gw_turbulence_integrals!(
    state::State,
    turbulence_scheme::Val{:TKEScheme},
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the gravity wave shear at ``(i, j, k)``.

The shear is given by 

```math 
S_\\mathrm{w} = \\sum_{r, \\lambda,\\mu,\\nu} m_r^2 \\left[\\frac{\\left|u_r\\right|^2 + \\left|v_r\\right|^2}{2}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu}
```

# Arguments

  - `state`: Model state.

  - `turbulence_scheme`:  General turbulence parameterization configuration.

  - `fc`: Coriolis parameter.

  - `omir`: Gravity-wave intrinsic frequency.

  - `wnrk`: Zonal wavenumber.

  - `wnrl`: Meridional wavenumber.

  - `wnrm`: Vertical wavenumber.

  - `wadr`: Phase-space wave-action density.

  - `n2r`: Squared buoyancy frequency at the ray volume position ``N_r^2``.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

"""
function compute_gw_turbulence_integrals! end

function compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    n2r::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; turbulence_scheme) = state.namelists.turbulence

    @dispatch_turbulence_scheme compute_gw_turbulence_integrals!(
        state,
        fc,
        omir,
        wnrk,
        wnrl,
        wnrm,
        wadr,
        factor,
        n2r,
        i,
        j,
        k,
        Val(turbulence_scheme),
    )
    return
end

function compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    n2r::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:NoTurbulence},
)
    return
end

@ivy function compute_gw_turbulence_integrals!(
    state::State,
    fc::AbstractFloat,
    omir::AbstractFloat,
    wnrk::AbstractFloat,
    wnrl::AbstractFloat,
    wnrm::AbstractFloat,
    wadr::AbstractFloat,
    factor::AbstractFloat,
    n2r::AbstractFloat,
    i::Integer,
    j::Integer,
    k::Integer,
    turbulence_scheme::Val{:TKEScheme},
)
    (; gwshear, gwbuoy) = state.turbulence.turbulencewkbintegrals
    (; wave_impact) = state.namelists.turbulence
    (; rhobar) = state.atmosphere

    if !wave_impact
        return
    end

    gwshear[i, j, k] +=
        wnrm^4 / omir * (fc^2 + omir^2) / (wnrk^2 + wnrl^2 + wnrm^2) * wadr /
        rhobar[i, j, k]

    bhat = sqrt(
        n2r^2 * (wnrk^2 + wnrl^2) / (wnrk^2 + wnrl^2 + wnrm^2) *
        2 *
        abs(wadr / factor / omir) / rhobar[i, j, k],
    )

    dphi = 2 * pi / 20
    phi = 0.0
    integral = 0.0
    while phi <= 2 * pi
        integral +=
            max(0, -(n2r + real(1im * wnrm * bhat * exp(1im * phi)))) * dphi
        phi += dphi
    end
    integral /= (2 * pi)
    gwbuoy[i, j, k] += integral * factor

    return
end

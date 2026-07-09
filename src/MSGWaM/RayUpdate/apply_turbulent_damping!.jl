"""
```julia
apply_turbulent_damping!(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    zr::AbstractFloat,
    dt::AbstractFloat,
)
```

Damping of the wave-action density due to turbulence.

The update of the physical-space wave-action density reads

```math
\\mathcal{A}_r \\rightarrow \\left[1 - 2\\Delta t \\left(\\gamma_s + \\gamma_w + \\gamma_w' \\right)\\right]\\mathcal{A}_r,
```

where the turbulent damping terms are given by

```math
\\begin{align*}
    \\gamma_s & = m_r^2 \\left[l_v\\left(1-\\delta_r\\right) + l_b\\delta_r\\right]\\Re\\left(Q_{0,r}\\right),\\\\
    \\gamma_w & = \\frac{m_r^2}{4} \\frac{N_r^2\\left(k_r^2+l_r^2\\right)}{N_r^2\\left(k_r^2+l_r^2\\right)+f^2m_r^2}\\left[l_v \\left(1-\\frac{f^2}{N_r^2}\\right)\\left(1+\\frac{k_r^2+l_r^2}{m_r^2}\\right)^{-1}-l_b\\right]\\Re\\left(Q_{2,r}\\right),\\\\
    \\gamma_w' & = -l_b\\frac{m_r}{2\\hat{\\omega}_r}\\sqrt{\\frac{N_r^2\\left(k_r^2+l_r^2\\right)}{\\left|\\boldsymbol{k}_r\\right|^2}\\frac{\\bar{\\rho}\\hat{\\omega}_r}{2\\mathcal{A}_r}}\\Re\\left(iQ_{1,r}\\right),
\\end{align*}
```

with

```math
\\delta_r = \\frac{N_r^2\\left(k_r^2+l_r^2\\right)}{2\\left[N_r^2\\left(k_r^2+l_r^2\\right)+f^2m_r^2\\right]}
```

and the turbulent mixing lengths ``l_v`` and ``l_b`` stored in `state.turbulence.turbulenceconstants.lv` and `state.turbulence.turbulenceconstants.lb`, respectively. Furthermore, the characteristic turbulent velocities ``Q_{0,r}``, ``Q_{1,r}`` and ``Q_{2,r}`` are computed with `compute_turbulent_velocity`.

# Arguments

  - `state`: Model state.

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `zr`: Position of the ray volume in ``z``.

  - `dt`: Time step.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_turbulent_velocity`](@ref)

!!! danger "Experimental"
    The turbulent damping of wave-action density is an experimental feature that hasn't been validated yet.
"""
function apply_turbulent_damping! end

@ivy function apply_turbulent_damping!(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    zr::AbstractFloat,
    dt::AbstractFloat,
)
    (; rays) = state.wkb
    (; lv, lb) = state.turbulence.turbulenceconstants
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; x_size, y_size) = state.namelists.domain
    (; rhobar) = state.atmosphere
    (; turbulent_damping) = state.namelists.wkb
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands

    if !turbulent_damping
        return
    end
    if rays.dens[r, i, j, k] == 0.0
        return
    end

    fc = coriolis_frequency * tref

    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]

    dkr = rays.dkray[r, i, j, k]
    dlr = rays.dlray[r, i, j, k]
    dmr = rays.dmray[r, i, j, k]

    kh2 = kr^2 + lr^2

    n2r = interpolate_stratification(zr, state, N2())

    omir = -sqrt(n2r * kh2 + fc^2 * mr^2) / sqrt(kh2 + mr^2)
    rhob = rhobar[i, j, k]

    factor = dmr

    if x_size > 1
        factor *= dkr
    end
    if y_size > 1
        factor *= dlr
    end

    wadr = rays.dens[r, i, j, k] * factor

    q1r = compute_turbulent_velocity(state, r, i, j, k, 1.0)
    q2r = compute_turbulent_velocity(state, r, i, j, k, 2.0)

    delta = n2r * kh2 / (2 * (n2r * kh2 + fc^2 * mr^2))

    gammas = 
        mr^2 * sqrt(2 * tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])) * 
        (lv * (1 - delta) + lb * delta)

    gammaw =
        mr^2 / 4 * n2r * kh2 / (n2r * kh2 + fc^2 * mr^2) *
        (lv * (1 - fc^2 / n2r) / (1 + kh2 / mr^2) - lb) *
        real(q2r)

    gammawp =
        -lb * mr / omir / 2 *
        sqrt(n2r^2 * kh2 / (kh2 + mr^2) * rhob / 2 / abs(omir / wadr)) *
        real(1im * q1r)

    wadr *= 1 - 2 * dt * (gammas + gammaw + gammawp)

    rays.dens[r, i, j, k] = wadr / factor

    return
end

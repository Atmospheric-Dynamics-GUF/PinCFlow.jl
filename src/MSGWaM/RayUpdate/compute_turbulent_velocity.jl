"""
```julia 
compute_turbulent_velocity(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::Tuple{<:Complex, <:Complex, <:Complex}
```

Compute and return the characteristic mean turbulent velocity amplitudes ``Q_{0,r}``, ``Q_{1,r}``, and ``Q_{2,r}`` for each ray volume.

```julia 
compute_turbulent_velocity(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    beta::AbstractFloat,
)::Complex
```

Compute and return the characteristic mean turbulent velocity amplitude ``Q_{\\beta,r}``, with ``\\beta`` given by the input parameter `beta`.

The velocity amplitude is computed from the numerical phase average 

```math 
Q_{0,r} = \\frac{1}{2\\pi}\\sum_{n=0}^{N_{\\phi}}\\sqrt{\\tilde{Q}_r^2(n\\Delta\\phi)}\\Delta\\phi\\;,
```

and for ``\\beta>0``

```math 
Q_{\\beta,r} = \\frac{1}{\\pi}\\sum_{n=0}^{N_{\\phi}}\\sqrt{\\tilde{Q}_r^2(n\\Delta\\phi)}e^{-i\\beta n\\Delta\\phi}\\Delta\\phi\\;.
```

```julia 
compute_turbulent_velocity(
    state::State,
    rhob::AbstractFloat,
    wadr::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    omir::AbstractFloat,
    phi::AbstractFloat,
)::AbstractFloat
```

Compute and return ``\\tilde{Q}_r`` with ``\\tilde{Q}_r^2`` being the leading-order turbulence contribution by shear production and buoyancy term given by

```math 
\\tilde{Q}_r^2 = \\max\\left\\{0,l_d\\left\\{l_v \\frac{m_r^2}{2}\\left[\\left|\\hat{\\mathbf{u}}_r\\right|^2-\\real\\left(\\hat{\\mathbf{u}}_r\\cdot\\hat{\\mathbf{u}}_r\\exp i2\\phi \\right)\\right]
-l_b\\left[N_r^2+\\real\\left(im_r\\hat{b}_r\\exp i\\phi\\right)\\right]\\right\\}\\right\\}\\;,
```

with

```math 
\\begin{align*}
\\left|\\hat{\\mathbf{u}}_r\\right|^2 &= \\frac{m_r^2 \\left(\\hat{\\omega}_r^2-f^2\\right)}{k_r^2+l_r^2+m_r^2}\\frac{2\\mathcal{A}_r}{\\hat{\\omega}_r\\bar{\\rho}} \\;, \\\\
\\hat{\\mathbf{u}}_r\\cdot\\hat{\\mathbf{u}}_r &= -\\frac{\\left(N_r^2+f^2\\right)\\left(k_r^2+l_r^2\\right)m_r^2}{\\left(k_r^2+l_r^2+m_r^2\\right)^2}\\frac{2\\mathcal{A}_r}{\\hat{\\omega}_r\\bar{\\rho}} \\;, \\\\
\\hat{b}_r &= \\sqrt{\\frac{N_r^2\\left(k_r^2+l_r^2\\right)}{\\left(k_r^2+l_r^2+m_r^2\\right)}\\frac{2\\mathcal{A}_r}{\\hat{\\omega}_r\\bar{\\rho}}} \\;,
\\end{align*}
```

and turbulence mixing lengths ``l_d``, ``l_v``, and ``l_b`` stored in `state.turbulence.turbulenceconstants.ld`, `state.turbulence.turbulenceconstants.lv`, and `state.turbulence.turbulenceconstants.lb`, respectively.

# Arguments 

  - `state`: Model state.

  - `r`: Ray-volume index. 

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `beta`: Index ``\\beta`` of ``Q_{\\beta,r}``.

  - `rhob`: Background density ``\\bar{\\rho}`` located at cell index ``(i,j,k)``.

  - `wadr`: Physical-space wave-action density ``\\mathcal{A}_r``. 

  - `kr`: Zonal wavenumber ``k_r``.

  - `lr`: Meridional wavenumber ``l_r``.

  - `mr`: Vertical wavenumber ``m_r``.

  - `n2r`: Buoyancy frequency at the ray volume position ``N_r^2``.

  - `fc`: Coriolis frequency ``f``.

  - `omir`: Intrinsic frequency ``\\hat{\\omega}_r``.

  - `phi`: Gravity wave phase ``\\phi``.
"""
function compute_turbulent_velocity end

function compute_turbulent_velocity(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::Tuple{<:Complex, <:Complex, <:Complex}
    (; rays) = state.wkb

    if rays.dens[r, i, j, k] == 0.0
        return 0.0, 0.0, 0.0
    end

    q00 = compute_turbulent_velocity(state, r, i, j, k, 0.0)
    q10 = compute_turbulent_velocity(state, r, i, j, k, 1.0)
    q20 = compute_turbulent_velocity(state, r, i, j, k, 2.0)

    return q00, q10, q20
end

function compute_turbulent_velocity(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    beta::AbstractFloat,
)::Complex
    (; rays) = state.wkb
    (; coriolis_frequency) = state.namelists.atmosphere
    (; tref) = state.constants
    (; rhobar) = state.atmosphere
    (; x_size, y_size) = state.namelists.domain
    (; branch) = state.namelists.wkb

    (xr, yr, zr) = get_physical_position(rays, r, i, j, k)

    rhob = rhobar[i, j, k]
    kr = rays.k[r, i, j, k]
    lr = rays.l[r, i, j, k]
    mr = rays.m[r, i, j, k]
    dkr = rays.dkray[r, i, j, k]
    dlr = rays.dlray[r, i, j, k]
    dmr = rays.dmray[r, i, j, k]
    n2r = interpolate_stratification(zr, state, N2())
    fc = coriolis_frequency * tref

    khr = sqrt(kr^2 + lr^2)

    omir = branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

    wadr = rays.dens[r, i, j, k] * dmr

    if x_size > 1
        wadr *= dkr
    end
    if y_size > 1
        wadr *= dlr
    end

    int_min = 0.0
    int_max = 2 * pi
    nphi = 20
    dphi = (int_max - int_min) / nphi
    phi = 0.0
    integral = 0.0
    while phi <= int_max
        qtilde = compute_turbulent_velocity(state, rhob, wadr, kr, lr, mr, n2r, fc, omir, phi)
        integral += qtilde * exp(-1im * beta * phi) * dphi
        phi += dphi
    end
    if beta == 0.0
        return integral / (2 * pi)
    else
        return integral / pi
    end
end

function compute_turbulent_velocity(
    state::State,
    rhob::AbstractFloat,
    wadr::AbstractFloat,
    kr::AbstractFloat,
    lr::AbstractFloat,
    mr::AbstractFloat,
    n2r::AbstractFloat,
    fc::AbstractFloat,
    omir::AbstractFloat,
    phi::AbstractFloat,
)::AbstractFloat
    (; ld, lv, lb) = state.turbulence.turbulenceconstants

    uhat2 =
        2 * mr^2 * wadr * (omir^2 + fc^2) / rhob / omir / (kr^2 + lr^2 + mr^2)
    u01u01 =
        -(n2r - fc^2) * (kr^2 + lr^2) * mr^2 / (kr^2 + lr^2 + mr^2)^2 *
        2 *
        wadr / omir / rhob
    bhat = sqrt(
        n2r^2 * (kr^2 + lr^2) / (kr^2 + lr^2 + mr^2) * 2 * abs(wadr / omir) / rhob,
    )

    sterm = mr^2 / 2 * (uhat2 - real(u01u01 * exp(2im * phi)))
    bterm = n2r + real(1im * mr * bhat * exp(1im * phi))

    qtilde = sqrt(max(0, ld * (lv * sterm - lb * bterm)))

    return qtilde
end

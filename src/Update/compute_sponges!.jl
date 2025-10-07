"""
```julia
compute_sponges!(state::State, dt::AbstractFloat)
```

Compute the Rayleigh-damping coefficients of the two sponges.

This method directly computes the Rayleigh-damping coefficient

```math
\\beta_\\mathrm{R} \\left(z\\right) = \\begin{cases}
    \\frac{\\beta_{\\mathrm{R}, \\max}}{\\Delta t} \\sin^2 \\left[\\frac{\\pi \\left(z - z_\\mathrm{R}\\right)}{2 \\left(L_z - z_\\mathrm{R}\\right)}\\right] & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
    0 & \\mathrm{else},
\\end{cases}
```

where ``\\beta_{\\mathrm{R}, \\max}`` and ``z_\\mathrm{R}`` are given by `state.namelists.sponge.betarmax` and `state.sponge.zsponge`, respectively. This coefficient is only used in the prognostic equations for the horizontal wind (if `state.namelists.sponge.damp_horizontal_wind_on_rhs` is `true`) and the transformed vertical wind. The corresponding damping terms are integrated on the right-hand sides.

This method also dispatches to a specific method that computes the Rayleigh damping coefficient of the RHS sponge defined for `state.namelists.sponge.sponge_type`.

```julia
compute_sponges!(
    state::State,
    dt::AbstractFloat,
    sponge_type::ExponentialSponge,
)
```

Compute the Rayleigh-damping coefficient of an exponential sponge.

If `state.namelists.sponge.lateral_sponge` is `true`, the Rayleigh-damping coefficient is

```math
\\alpha_\\mathrm{R} \\left(x, y, z\\right) = \\frac{\\alpha_{\\mathrm{R}, \\max}}{3} \\left[\\exp \\left(\\frac{\\left|x\\right| - L_x / 2}{\\Delta x_\\mathrm{R}}\\right) + \\exp \\left(\\frac{\\left|y\\right| - L_y / 2}{\\Delta y_\\mathrm{R}}\\right) + \\exp \\left(\\frac{z - L_z}{\\Delta z_\\mathrm{R}}\\right)\\right]
```

and otherwise, it is

```math
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, \\max} \\left[\\exp \\left(\\frac{z - L_z}{\\Delta z_\\mathrm{R}}\\right)\\right],
```

where ``\\alpha_{\\mathrm{R}, \\max}``, ``\\Delta x_\\mathrm{R}``, ``\\Delta y_\\mathrm{R}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.alpharmax`, `state.sponge.dxsponge`, `state.sponge.dysponge` and `state.sponge.dzsponge`, respectively. If the grid size in a horizontal direction is one, the contribution from that direction is set to zero and the other two are reweighted accordingly.

```julia
compute_sponges!(state::State, dt::AbstractFloat, sponge_type::COSMOSponge)
```

Compute the Rayleigh-damping coefficient of a sponge similar to that used by the COSMO model.

If `state.namelists.sponge.lateral_sponge` is `true`, the Rayleigh-damping coefficient is

```math
\\alpha_\\mathrm{R} \\left(x, y, z\\right) = \\alpha_{\\mathrm{R}, x} \\left(x\\right) + \\alpha_{\\mathrm{R}, y} \\left(y\\right) + \\alpha_{\\mathrm{R}, z} \\left(z\\right)
```

with

```math
\\begin{align*}
    \\alpha_{\\mathrm{R}, x} \\left(x\\right) & = \\begin{cases}
        \\left(2 N_\\mathrm{R} \\Delta t\\right)^{- 1} \\left\\{1 - \\cos \\left[\\frac{\\pi \\left(x_{\\mathrm{R}, 0} - x\\right)}{\\Delta x_\\mathrm{R}}\\right]\\right\\} & \\mathrm{if} \\quad x \\leq x_{\\mathrm{R}, 0},\\\\
        \\left(2 N_\\mathrm{R} \\Delta t\\right)^{- 1} \\left\\{1 - \\cos \\left[\\frac{\\pi \\left(x - x_{\\mathrm{R}, 1}\\right)}{\\Delta x_\\mathrm{R}}\\right]\\right\\} & \\mathrm{if} \\quad x \\geq x_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, y} \\left(y\\right) & = \\begin{cases}
        \\left(2 N_\\mathrm{R} \\Delta t\\right)^{- 1} \\left\\{1 - \\cos \\left[\\frac{\\pi \\left(y_{\\mathrm{R}, 0} - y\\right)}{\\Delta y_\\mathrm{R}}\\right]\\right\\} & \\mathrm{if} \\quad y \\leq y_{\\mathrm{R}, 0},\\\\
        \\left(2 N_\\mathrm{R} \\Delta t\\right)^{- 1} \\left\\{1 - \\cos \\left[\\frac{\\pi \\left(y - y_{\\mathrm{R}, 1}\\right)}{\\Delta y_\\mathrm{R}}\\right]\\right\\} & \\mathrm{if} \\quad y \\geq y_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, z} \\left(z\\right) & = \\begin{cases}
        \\left(2 N_\\mathrm{R} \\Delta t\\right)^{- 1} \\left\\{1 - \\cos \\left[\\frac{\\pi \\left(z - z_\\mathrm{R}\\right)}{\\Delta z_\\mathrm{R}}\\right]\\right\\} & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
        0 & \\mathrm{else}
    \\end{cases}
\\end{align*}
```

and otherwise, it is

```math
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right),
```

where ``N_\\mathrm{R}``, ``x_{\\mathrm{R}, 0}``, ``x_{\\mathrm{R}, 1}``, ``y_{\\mathrm{R}, 0}``, ``y_{\\mathrm{R}, 1}``, ``z_{\\mathrm{R}, 0}`` and ``z_{\\mathrm{R}, 1}`` are given by `state.namelists.sponge.cosmo_steps` and the properties `xsponge0`, `xsponge1`, `ysponge0`, `ysponge1`, `zsponge0` and `zsponge1` of `state.sponge`, respectively.

```julia
compute_sponges!(state::State, dt::AbstractFloat, sponge_type::PolynomialSponge)
```

Compute the Rayleigh-damping coefficient of a polynomial sponge.

If `state.namelists.sponge.lateral_sponge` is `true`, the Rayleigh-damping coefficient is

```math
\\alpha_\\mathrm{R} \\left(x, y, z\\right) = \\frac{\\alpha_{\\mathrm{R}, x} \\left(x\\right) + \\alpha_{\\mathrm{R}, y} \\left(y\\right) + \\alpha_{\\mathrm{R}, z} \\left(z\\right)}{3}
```

with

```math
\\begin{align*}
    \\alpha_{\\mathrm{R}, x} \\left(x\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\left(\\frac{x_{\\mathrm{R}, 0} - x}{\\Delta x_\\mathrm{R}}\\right)^{n_\\mathrm{R}} & \\mathrm{if} \\quad x \\leq x_{\\mathrm{R}, 0},\\\\
        \\alpha_{\\mathrm{R}, \\max} \\left(\\frac{x - x_{\\mathrm{R}, 1}}{\\Delta x_\\mathrm{R}}\\right)^{n_\\mathrm{R}} & \\mathrm{if} \\quad x \\geq x_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, y} \\left(y\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\left(\\frac{y_{\\mathrm{R}, 0} - y}{\\Delta y_\\mathrm{R}}\\right)^{n_\\mathrm{R}} & \\mathrm{if} \\quad y \\leq y_{\\mathrm{R}, 0},\\\\
        \\alpha_{\\mathrm{R}, \\max} \\left(\\frac{y - y_{\\mathrm{R}, 1}}{\\Delta y_\\mathrm{R}}\\right)^{n_\\mathrm{R}} & \\mathrm{if} \\quad y \\geq y_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, z} \\left(z\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\left(\\frac{z - z_\\mathrm{R}}{\\Delta z_\\mathrm{R}}\\right)^{n_\\mathrm{R}} & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
        0 & \\mathrm{else}
    \\end{cases}
\\end{align*}
```

and otherwise, it is

```math
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right),
```

where ``n_\\mathrm{R}`` is given by `state.namelists.sponge.sponge_order`. If the grid size in a horizontal direction is one, the contribution from that direction is set to zero and the other two are reweighted accordingly.

```julia
compute_sponges!(state::State, dt::AbstractFloat, sponge_type::SinusoidalSponge)
```

Compute the Rayleigh-damping coefficient of a sinusoidal sponge.

If `state.namelists.sponge.lateral_sponge` is `true`, the Rayleigh-damping coefficient is

```math
\\alpha_\\mathrm{R} \\left(x, y, z\\right) = \\frac{\\alpha_{\\mathrm{R}, x} \\left(x\\right) + \\alpha_{\\mathrm{R}, y} \\left(y\\right) + \\alpha_{\\mathrm{R}, z} \\left(z\\right)}{3}
```

with

```math
\\begin{align*}
    \\alpha_{\\mathrm{R}, x} \\left(x\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\sin^2 \\left[\\frac{\\pi \\left(x_{\\mathrm{R}, 0} - x\\right)}{2 \\Delta x_\\mathrm{R}}\\right] & \\mathrm{if} \\quad x \\leq x_{\\mathrm{R}, 0},\\\\
        \\alpha_{\\mathrm{R}, \\max} \\sin^2 \\left[\\frac{\\pi \\left(x - x_{\\mathrm{R}, 1}\\right)}{2 \\Delta x_\\mathrm{R}}\\right] & \\mathrm{if} \\quad x \\geq x_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, y} \\left(y\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\sin^2 \\left[\\frac{\\pi \\left(y_{\\mathrm{R}, 0} - y\\right)}{2 \\Delta y_\\mathrm{R}}\\right] & \\mathrm{if} \\quad y \\leq y_{\\mathrm{R}, 0},\\\\
        \\alpha_{\\mathrm{R}, \\max} \\sin^2 \\left[\\frac{\\pi \\left(y - y_{\\mathrm{R}, 1}\\right)}{2 \\Delta y_\\mathrm{R}}\\right] & \\mathrm{if} \\quad y \\geq y_{\\mathrm{R}, 1},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_{\\mathrm{R}, z} \\left(z\\right) & = \\begin{cases}
        \\alpha_{\\mathrm{R}, \\max} \\sin^2 \\left[\\frac{\\pi \\left(z - z_\\mathrm{R}\\right)}{2 \\Delta z_\\mathrm{R}}\\right] & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
        0 & \\mathrm{else}
    \\end{cases}
\\end{align*}
```

and otherwise, it is

```math
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right).
```

If the grid size in a horizontal direction is one, the contribution from that direction is set to zero and the other two are reweighted accordingly.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `sponge_type`: Specification of the spatial dependence of the  Rayleigh-damping coefficient.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_sponges! end

function compute_sponges!(state::State, dt::AbstractFloat)
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; zc, lz) = state.grid
    (; betar, zsponge) = state.sponge
    (; use_sponge, sponge_type, betarmax) = state.namelists.sponge

    if !use_sponge
        return
    end

    compute_sponges!(state, dt, sponge_type)

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        if zc[i, j, k] >= zsponge
            betar[i, j, k] =
                betarmax / dt *
                sin(0.5 * pi * (zc[i, j, k] - zsponge) / (lz - zsponge))^2.0
        end
    end

    @ivy ko == 0 && betar[:, :, k0 - 1] .= betar[:, :, k0]
    @ivy ko + nzz == zz_size && betar[:, :, k1 + 1] .= betar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    sponge_type::ExponentialSponge,
)
    (; namelists, domain) = state
    (; x_size, y_size, z_size) = namelists.domain
    (; zz_size, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, zc, lx, ly, lz) = state.grid
    (; tref) = state.constants
    (; lateral_sponge, alpharmax) = namelists.sponge
    (; alphar, dxsponge, dysponge, dzsponge) = state.sponge

    alpharzmax = alpharmax * tref

    dim = 1
    x_size > 1 && (dim += 1)
    y_size > 1 && (dim += 1)
    lateral_sponge && (alpharzmax /= dim)
    alpharxmax = alpharymax = alpharzmax

    alphar .= 0.0

    imin = lateral_sponge ? i0 : 1
    imax = lateral_sponge ? i1 : nxx

    jmin = lateral_sponge ? j0 : 1
    jmax = lateral_sponge ? j1 : nyy

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in jmin:jmax, i in imin:imax
        height = zc[i, j, k]

        if z_size > 1
            alphar[i, j, k] =
                alphar[i, j, k] + alpharzmax * exp((height - lz) / dzsponge)
        end
        if lateral_sponge
            if x_size > 1
                if x[io + i] <= 0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax * exp((-lx / 2 - x[io + i]) / dxsponge)
                else
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax * exp((x[io + i] - lx / 2) / dxsponge)
                end
            end
            if y_size > 1
                if y[jo + j] <= 0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax * exp((-ly / 2 - y[jo + j]) / dysponge)
                else
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax * exp((y[jo + j] - ly / 2) / dysponge)
                end
            end
        end
    end

    if lateral_sponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    @ivy ko == 0 && alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    @ivy ko + nzz == zz_size && alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    sponge_type::COSMOSponge,
)
    (; namelists, domain) = state
    (; x_size, y_size, z_size) = namelists.domain
    (; zz_size, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, zc) = state.grid
    (; lateral_sponge, cosmo_steps) = namelists.sponge
    (;
        alphar,
        dxsponge,
        dysponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        zsponge,
    ) = state.sponge

    alphar .= 0.0

    imin = lateral_sponge ? i0 : 1
    imax = lateral_sponge ? i1 : nxx

    jmin = lateral_sponge ? j0 : 1
    jmax = lateral_sponge ? j1 : nyy

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in jmin:jmax, i in imin:imax
        height = zc[i, j, k]

        if z_size > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    0.5 / cosmo_steps / dt *
                    (1.0 - cos(pi * (height - zsponge) / dzsponge))
            end
        end
        if lateral_sponge
            if x_size > 1
                if x[io + i] <= xsponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmo_steps / dt *
                        (1.0 - cos(pi * (xsponge0 - x[io + i]) / dxsponge))
                elseif x[io + i] >= xsponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmo_steps / dt *
                        (1.0 - cos(pi * (x[io + i] - xsponge1) / dxsponge))
                end
            end
            if y_size > 1
                if y[jo + j] <= ysponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmo_steps / dt *
                        (1.0 - cos(pi * (ysponge0 - y[jo + j]) / dysponge))
                elseif y[jo + j] >= ysponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmo_steps / dt *
                        (1.0 - cos(pi * (y[jo + j] - ysponge1) / dysponge))
                end
            end
        end
    end

    if lateral_sponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    @ivy ko == 0 && alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    @ivy ko + nzz == zz_size && alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    sponge_type::PolynomialSponge,
)
    (; namelists, domain) = state
    (; x_size, y_size, z_size) = namelists.domain
    (; zz_size, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, zc) = state.grid
    (; tref) = state.constants
    (; lateral_sponge, alpharmax, sponge_order) = namelists.sponge
    (;
        alphar,
        dxsponge,
        dysponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        zsponge,
    ) = state.sponge

    alpharzmax = alpharmax * tref

    dim = 1
    x_size > 1 && (dim += 1)
    y_size > 1 && (dim += 1)
    lateral_sponge && (alpharzmax /= dim)
    alpharxmax = alpharymax = alpharzmax

    alphar .= 0.0

    imin = lateral_sponge ? i0 : 1
    imax = lateral_sponge ? i1 : nxx

    jmin = lateral_sponge ? j0 : 1
    jmax = lateral_sponge ? j1 : nyy

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in jmin:jmax, i in imin:imax
        height = zc[i, j, k]

        if z_size > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    alpharzmax * ((height - zsponge) / dzsponge)^sponge_order
            end
        end
        if lateral_sponge
            if x_size > 1
                if x[io + i] <= xsponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        ((xsponge0 - x[io + i]) / dxsponge)^sponge_order
                elseif x[io + i] >= xsponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        ((x[io + i] - xsponge1) / dxsponge)^sponge_order
                end
            end
            if y_size > 1
                if y[jo + j] <= ysponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        ((ysponge0 - y[jo + j]) / dysponge)^sponge_order
                elseif y[jo + j] >= ysponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        ((y[jo + j] - ysponge1) / dysponge)^sponge_order
                end
            end
        end
    end

    if lateral_sponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    @ivy ko == 0 && alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    @ivy ko + nzz == zz_size && alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    sponge_type::SinusoidalSponge,
)
    (; namelists, domain) = state
    (; x_size, y_size, z_size) = namelists.domain
    (; zz_size, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, zc) = state.grid
    (; tref) = state.constants
    (; lateral_sponge, alpharmax) = namelists.sponge
    (;
        alphar,
        dxsponge,
        dysponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        zsponge,
    ) = state.sponge

    alpharzmax = alpharmax * tref

    dim = 1
    x_size > 1 && (dim += 1)
    y_size > 1 && (dim += 1)
    lateral_sponge && (alpharzmax /= dim)
    alpharxmax = alpharymax = alpharzmax

    alphar .= 0.0

    imin = lateral_sponge ? i0 : 1
    imax = lateral_sponge ? i1 : nxx

    jmin = lateral_sponge ? j0 : 1
    jmax = lateral_sponge ? j1 : nyy

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in jmin:jmax, i in imin:imax
        height = zc[i, j, k]

        if z_size > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    alpharzmax *
                    sin(0.5 * pi * (height - zsponge) / dzsponge)^2.0
            end
        end
        if lateral_sponge
            if x_size > 1
                if x[io + i] <= xsponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        sin(0.5 * pi * (xsponge0 - x[io + i]) / dxsponge)^2.0
                elseif x[io + i] >= xsponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        sin(0.5 * pi * (x[io + i] - xsponge1) / dxsponge)^2.0
                end
            end
            if y_size > 1
                if y[jo + j] <= ysponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        sin(0.5 * pi * (ysponge0 - y[jo + j]) / dysponge)^2.0
                elseif y[jo + j] >= ysponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        sin(0.5 * pi * (y[jo + j] - ysponge1) / dysponge)^2.0
                end
            end
        end
    end

    if lateral_sponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    @ivy ko == 0 && alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    @ivy ko + nzz == zz_size && alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

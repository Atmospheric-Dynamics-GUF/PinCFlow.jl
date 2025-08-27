"""
```julia
compute_sponges!(state::State, dt::AbstractFloat)
```

Compute the Rayleigh-damping coefficients of the two sponges.

This method directly computes the Rayleigh-damping coefficients

```math
\\beta_\\mathrm{R}^{uv} = \\beta_\\mathrm{R}^{\\widehat{w}} \\left(z\\right) = \\begin{cases}
    \\frac{\\beta_{\\mathrm{R}, \\max}}{\\Delta t} \\sin^2 \\left[\\frac{\\pi \\left(z - z_\\mathrm{R}\\right)}{2 \\left(L_z - z_\\mathrm{R}\\right)}\\right] & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
    0 & \\mathrm{else},
\\end{cases}
```

where ``\\beta_{\\mathrm{R}, \\max}`` and ``z_\\mathrm{R}`` are given by `state.namelists.sponge.betarmax` and `state.sponge.zsponge`, respectively. These coefficients are only used in the prognostic equations for the horizontal wind (``\\beta_\\mathrm{R}^{uv}``, only if `state.namelists.sponge.sponge_uv` is `true`) and the transformed vertical wind (``\\beta_\\mathrm{R}^{\\widehat{w}}``). The corresponding damping terms are integrated on the right-hand sides.

This method also dispatches to a specific method that computes the Rayleigh damping coefficient of the RHS sponge defined for `state.namelists.sponge.spongetype`.

```julia
compute_sponges!(state::State, dt::AbstractFloat, spongetype::ExponentialSponge)
```

Compute the Rayleigh-damping coefficient of an exponential sponge.

If `state.namelists.sponge.lateralsponge` is `true`, the Rayleigh-damping coefficient is

```math
\\alpha_\\mathrm{R} \\left(x, y, z\\right) = \\frac{\\alpha_{\\mathrm{R}, \\max}}{3} \\left[\\exp \\left(\\frac{\\left|x\\right| - L_x / 2}{\\Delta x_\\mathrm{R}}\\right) + \\exp \\left(\\frac{\\left|y\\right| - L_y / 2}{\\Delta y_\\mathrm{R}}\\right) + \\exp \\left(\\frac{z - L_z}{\\Delta z_\\mathrm{R}}\\right)\\right]
```

and otherwise, it is

```math
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, \\max} \\left[\\exp \\left(\\frac{z - L_z}{\\Delta z_\\mathrm{R}}\\right)\\right],
```

where ``\\alpha_{\\mathrm{R}, \\max}``, ``\\Delta x_\\mathrm{R}``, ``\\Delta y_\\mathrm{R}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.alpharmax`, `state.sponge.dxsponge`, `state.sponge.dysponge` and `state.sponge.dzsponge`, respectively. If the grid size in a horizontal direction is one, the contribution from that direction is set to zero and the other two are reweighted accordingly.

```julia
compute_sponges!(state::State, dt::AbstractFloat, spongetype::COSMOSponge)
```

Compute the Rayleigh-damping coefficient of a sponge similar to that used by the COSMO model.

If `state.namelists.sponge.lateralsponge` is `true`, the Rayleigh-damping coefficient is

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

where ``N_\\mathrm{R}``, ``x_{\\mathrm{R}, 0}``, ``x_{\\mathrm{R}, 1}``, ``y_{\\mathrm{R}, 0}``, ``y_{\\mathrm{R}, 1}``, ``z_{\\mathrm{R}, 0}`` and ``z_{\\mathrm{R}, 1}`` are given by `state.namelists.sponge.cosmosteps` and the properties `xsponge0`, `xsponge1`, `ysponge0`, `ysponge1`, `zsponge0` and `zsponge1` of `state.sponge`, respectively.

```julia
compute_sponges!(state::State, dt::AbstractFloat, spongetype::PolynomialSponge)
```

Compute the Rayleigh-damping coefficient of a polynomial sponge.

If `state.namelists.sponge.lateralsponge` is `true`, the Rayleigh-damping coefficient is

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

where ``n_\\mathrm{R}`` is given by `state.namelists.sponge.spongeorder`. If the grid size in a horizontal direction is one, the contribution from that direction is set to zero and the other two are reweighted accordingly.

```julia
compute_sponges!(state::State, dt::AbstractFloat, spongetype::SinusoidalSponge)
```

Compute the Rayleigh-damping coefficient of a sinusoidal sponge.

If `state.namelists.sponge.lateralsponge` is `true`, the Rayleigh-damping coefficient is

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

  - `spongetype`: Specification of the spatial dependence of the  Rayleigh-damping coefficient.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_sponges! end

function compute_sponges!(state::State, dt::AbstractFloat)
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; ztfc, lz, jac) = state.grid
    (; betar, zsponge) = state.sponge
    (; spongelayer, spongetype, betarmax) = state.namelists.sponge

    if !spongelayer
        return
    end

    compute_sponges!(state, dt, spongetype)

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
        if ztfc[i, j, k] >= zsponge
            betar[i, j, k] =
                betarmax / dt *
                sin(0.5 * pi * (ztfc[i, j, k] - zsponge) / (lz - zsponge))^2.0
        end
    end

    if ko == 0
        @views betar[:, :, k0 - 1] .= betar[:, :, k0]
    end

    if ko + nzz == sizezz
        @views betar[:, :, k1 + 1] .= betar[:, :, k1]
    end

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    spongetype::ExponentialSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc, lx, ly, lz) = state.grid
    (; tref) = state.constants
    (; lateralsponge, alpharmax) = namelists.sponge
    (; alphar, dxsponge, dysponge, dzsponge) = state.sponge

    alpharzmax = alpharmax * tref

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            alpharzmax = alpharzmax / 3.0
            alpharxmax = alpharzmax
            alpharymax = alpharzmax
        elseif sizex > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = alpharzmax
            alpharymax = 0.0
        elseif sizey > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = 0.0
            alpharymax = alpharzmax
        else
            alpharxmax = 0.0
            alpharymax = 0.0
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphar .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            alphar[i, j, k] =
                alphar[i, j, k] + alpharzmax * exp((height - lz) / dzsponge)
        end
        if lateralsponge
            if sizex > 1
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
            if sizey > 1
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

    if lateralsponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    ko == 0 && @views alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    ko + nzz == sizezz && @views alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    spongetype::COSMOSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc) = state.grid
    (; lateralsponge, cosmosteps) = namelists.sponge
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

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphar .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    0.5 / cosmosteps / dt *
                    (1.0 - cos(pi * (height - zsponge) / dzsponge))
            end
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= xsponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (xsponge0 - x[io + i]) / dxsponge))
                elseif x[io + i] >= xsponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (x[io + i] - xsponge1) / dxsponge))
                end
            end
            if sizey > 1
                if y[jo + j] <= ysponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (ysponge0 - y[jo + j]) / dysponge))
                elseif y[jo + j] >= ysponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (y[jo + j] - ysponge1) / dysponge))
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    ko == 0 && @views alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    ko + nzz == sizezz && @views alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    spongetype::PolynomialSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc) = state.grid
    (; tref) = state.constants
    (; lateralsponge, alpharmax, spongeorder) = namelists.sponge
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

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            alpharzmax = alpharzmax / 3.0
            alpharxmax = alpharzmax
            alpharymax = alpharzmax
        elseif sizex > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = alpharzmax
            alpharymax = 0.0
        elseif sizey > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = 0.0
            alpharymax = alpharzmax
        else
            alpharxmax = 0.0
            alpharymax = 0.0
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphar .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    alpharzmax * ((height - zsponge) / dzsponge)^spongeorder
            end
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= xsponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        ((xsponge0 - x[io + i]) / dxsponge)^spongeorder
                elseif x[io + i] >= xsponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharxmax *
                        ((x[io + i] - xsponge1) / dxsponge)^spongeorder
                end
            end
            if sizey > 1
                if y[jo + j] <= ysponge0
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        ((ysponge0 - y[jo + j]) / dysponge)^spongeorder
                elseif y[jo + j] >= ysponge1
                    alphar[i, j, k] =
                        alphar[i, j, k] +
                        alpharymax *
                        ((y[jo + j] - ysponge1) / dysponge)^spongeorder
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    ko == 0 && @views alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    ko + nzz == sizezz && @views alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

function compute_sponges!(
    state::State,
    dt::AbstractFloat,
    spongetype::SinusoidalSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc) = state.grid
    (; tref) = state.constants
    (; lateralsponge, alpharmax) = namelists.sponge
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

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            alpharzmax = alpharzmax / 3.0
            alpharxmax = alpharzmax
            alpharymax = alpharzmax
        elseif sizex > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = alpharzmax
            alpharymax = 0.0
        elseif sizey > 1
            alpharzmax = alpharzmax / 2.0
            alpharxmax = 0.0
            alpharymax = alpharzmax
        else
            alpharxmax = 0.0
            alpharymax = 0.0
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphar .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphar[i, j, k] =
                    alphar[i, j, k] +
                    alpharzmax *
                    sin(0.5 * pi * (height - zsponge) / dzsponge)^2.0
            end
        end
        if lateralsponge
            if sizex > 1
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
            if sizey > 1
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

    if lateralsponge
        set_zonal_boundaries_of_field!(alphar, namelists, domain)
        set_meridional_boundaries_of_field!(alphar, namelists, domain)
    end

    ko == 0 && @views alphar[:, :, k0 - 1] .= alphar[:, :, k0]
    ko + nzz == sizezz && @views alphar[:, :, k1 + 1] .= alphar[:, :, k1]

    return
end

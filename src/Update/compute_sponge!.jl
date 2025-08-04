"""
```julia
compute_sponge!(state::State, dt::AbstractFloat)
```

Compute the Rayleigh-damping coefficient(s) of the sponge layer.

If `state.namelists.sponge.unifiedsponge` is `false`, this method directly computes the Rayleigh-damping coefficients

```math
\\begin{align*}
    \\alpha_\\mathrm{R}^{\\rho'} \\left(z\\right) & = \\begin{cases}
        \\frac{a_\\mathrm{R}}{\\Delta t} \\sin^2 \\left[\\frac{\\pi \\left(z - z_\\mathrm{R}\\right)}{2 \\left(L_z - z_\\mathrm{R}\\right)}\\right] & \\mathrm{if} \\quad z \\geq z_\\mathrm{R},\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    \\alpha_\\mathrm{R}^{\\widehat{w}} \\left(z\\right) & = \\frac{\\alpha_\\mathrm{R}^{\\rho'}}{J},
\\end{align*}
```

where ``a_\\mathrm{R}`` and ``z_\\mathrm{R}`` are given by `state.namelists.sponge.spongealphaz_fac` and `state.sponge.zsponge`, respectively. These coefficients are only used in the prognostic equations for the density fluctuations (``\\alpha_\\mathrm{R}^{\\rho'}``) and transformed vertical wind (``\\alpha_\\mathrm{R}^{\\widehat{w}}``). The corresponding damping terms are integrated on the right-hand sides. If `state.namelists.sponge.unifiedsponge` is `true`, this method dispatches to a specific method that computes the Rayleigh damping coefficient of the sponge defined for `state.namelists.sponge.spongetype`.

```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::ExponentialSponge)
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

where ``\\alpha_{\\mathrm{R}, \\max}``, ``\\Delta x_\\mathrm{R}``, ``\\Delta y_\\mathrm{R}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.spongealphaz_dim`, `state.sponge.dxsponge`, `state.sponge.dysponge` and `state.sponge.dzsponge`, respectively. In the notation used here, the horizontal boundaries of the domain are ``\\left(- L_x / 2, L_x / 2\\right)`` and ``\\left(- L_y / 2, L_y / 2\\right)``. In the general case where this is not necessarily true, the coefficient is defined such that its minimum is at the horizontal center of the domain.

```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::COSMOSponge)
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
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right)
```

where ``N_\\mathrm{R}``, ``x_{\\mathrm{R}, 0}``, ``x_{\\mathrm{R}, 1}``, ``\\Delta x_\\mathrm{R}``, ``y_{\\mathrm{R}, 0}``, ``y_{\\mathrm{R}, 1}``, ``\\Delta y_\\mathrm{R}``, ``z_{\\mathrm{R}, 0}``, ``z_{\\mathrm{R}, 1}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.cosmosteps`, and the properties `xsponge0`, `xsponge1`, `dxsponge`, `ysponge0`, `ysponge1`, `dysponge`, `zsponge0`, `zsponge1` and `dzsponge` of `state.sponge`, respectively.

```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::PolynomialSponge)
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
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right)
```

where ``\\alpha_{\\mathrm{R}, \\max}``, ``n_\\mathrm{R}``, ``x_{\\mathrm{R}, 0}``, ``x_{\\mathrm{R}, 1}``, ``\\Delta x_\\mathrm{R}``, ``y_{\\mathrm{R}, 0}``, ``y_{\\mathrm{R}, 1}``, ``\\Delta y_\\mathrm{R}``, ``z_{\\mathrm{R}, 0}``, ``z_{\\mathrm{R}, 1}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.spongealphaz_dim`, `state.namelists.sponge.spongeorder`, and the properties `xsponge0`, `xsponge1`, `dxsponge`, `ysponge0`, `ysponge1`, `dysponge`, `zsponge0`, `zsponge1` and `dzsponge` of `state.sponge`, respectively.

```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::SinusoidalSponge)
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
\\alpha_\\mathrm{R} \\left(z\\right) = \\alpha_{\\mathrm{R}, z} \\left(z\\right)
```

where ``\\alpha_{\\mathrm{R}, \\max}``, ``x_{\\mathrm{R}, 0}``, ``x_{\\mathrm{R}, 1}``, ``\\Delta x_\\mathrm{R}``, ``y_{\\mathrm{R}, 0}``, ``y_{\\mathrm{R}, 1}``, ``\\Delta y_\\mathrm{R}``, ``z_{\\mathrm{R}, 0}``, ``z_{\\mathrm{R}, 1}`` and ``\\Delta z_\\mathrm{R}`` are given by `state.namelists.sponge.spongealphaz_dim` and the properties `xsponge0`, `xsponge1`, `dxsponge`, `ysponge0`, `ysponge1`, `dysponge`, `zsponge0`, `zsponge1` and `dzsponge` of `state.sponge`, respectively.

# Arguments

  - `state`: Model state.
  - `dt`: Time step.
  - `spongetype`: Specification of the spatial dependence of the  Rayleigh-damping coefficient.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_sponge! end

function compute_sponge!(state::State, dt::AbstractFloat)
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; ztfc, lz, jac) = state.grid
    (; kr_sp_tfc, kr_sp_w_tfc, zsponge) = state.sponge
    (; unifiedsponge, spongetype, spongealphaz_fac) = state.namelists.sponge

    if unifiedsponge
        compute_sponge!(state, dt, spongetype)
    else
        alpspg = spongealphaz_fac / dt

        kz0 = ko == 0 ? k0 : k0 - 1
        kz1 = ko + nzz == sizezz ? k1 : k1 + 1

        for k in kz0:kz1, j in (j0 - 1):(j1 + 1), i in (i0 - 1):(i1 + 1)
            if ztfc[i, j, k] >= zsponge
                kr_sp_tfc[i, j, k] =
                    alpspg *
                    sin(
                        0.5 * pi * (ztfc[i, j, k] - zsponge) /
                        (lz[2] - zsponge),
                    )^2.0
                kr_sp_w_tfc[i, j, k] = kr_sp_tfc[i, j, k] / jac[i, j, k]
            end
        end

        if ko == 0
            @views kr_sp_tfc[:, :, k0 - 1] .= kr_sp_tfc[:, :, k0]
            @views kr_sp_w_tfc[:, :, k0 - 1] .= kr_sp_w_tfc[:, :, k0]
        end

        if ko + nzz == sizezz
            @views kr_sp_tfc[:, :, k1 + 1] .= kr_sp_tfc[:, :, k1]
            @views kr_sp_w_tfc[:, :, k1 + 1] .= kr_sp_w_tfc[:, :, k1]
        end
    end

    return
end

function compute_sponge!(
    state::State,
    dt::AbstractFloat,
    spongetype::ExponentialSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc, lx, ly, lz) = state.grid
    (; tref) = state.constants
    (; lateralsponge, spongealphaz_dim) = namelists.sponge
    (; alphaunifiedsponge, dxsponge, dysponge, dzsponge) = state.sponge

    spongealphaz = spongealphaz_dim * tref

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            spongealphaz = spongealphaz / 3.0
            spongealphax = spongealphaz
            spongealphay = spongealphaz
        elseif sizex > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = spongealphaz
            spongealphay = 0.0
        elseif sizey > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = 0.0
            spongealphay = spongealphaz
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphaunifiedsponge .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            alphaunifiedsponge[i, j, k] =
                alphaunifiedsponge[i, j, k] +
                spongealphaz * exp((height - lz[2]) / dzsponge)
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= 0.5 * (lx[1] + lx[2])
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax * exp((lx[1] - x[io + i]) / dxsponge)
                else
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax * exp((x[io + i] - lx[2]) / dxsponge)
                end
            end
            if sizey > 1
                if y[jo + j] <= 0.5 * (ly[1] + ly[2])
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay * exp((ly[1] - y[jo + j]) / dysponge)
                else
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay * exp((y[jo + j] - ly[2]) / dysponge)
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphaunifiedsponge, namelists, domain)
        set_meridional_boundaries_of_field!(
            alphaunifiedsponge,
            namelists,
            domain,
        )
    end

    ko == 0 &&
        @views alphaunifiedsponge[:, :, k0 - 1] .= alphaunifiedsponge[:, :, k0]
    ko + nzz == sizezz &&
        @views alphaunifiedsponge[:, :, k1 + 1] .= alphaunifiedsponge[:, :, k1]

    return
end

function compute_sponge!(
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
        alphaunifiedsponge,
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

    alphaunifiedsponge .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphaunifiedsponge[i, j, k] =
                    alphaunifiedsponge[i, j, k] +
                    0.5 / cosmosteps / dt *
                    (1.0 - cos(pi * (height - zsponge) / dzsponge))
            end
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= xsponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (xsponge0 - x[io + i]) / dxsponge))
                elseif x[io + i] >= xsponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (x[io + i] - xsponge1) / dxsponge))
                end
            end
            if sizey > 1
                if y[jo + j] <= ysponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (ysponge0 - y[jo + j]) / dysponge))
                elseif y[jo + j] >= ysponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        0.5 / cosmosteps / dt *
                        (1.0 - cos(pi * (y[jo + j] - ysponge1) / dysponge))
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphaunifiedsponge, namelists, domain)
        set_meridional_boundaries_of_field!(
            alphaunifiedsponge,
            namelists,
            domain,
        )
    end

    ko == 0 &&
        @views alphaunifiedsponge[:, :, k0 - 1] .= alphaunifiedsponge[:, :, k0]
    ko + nzz == sizezz &&
        @views alphaunifiedsponge[:, :, k1 + 1] .= alphaunifiedsponge[:, :, k1]

    return
end

function compute_sponge!(
    state::State,
    dt::AbstractFloat,
    spongetype::PolynomialSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc) = state.grid
    (; tref) = state.constants
    (; lateralsponge, spongealphaz_dim, spongeorder) = namelists.sponge
    (;
        alphaunifiedsponge,
        dxsponge,
        dysponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        zsponge,
    ) = state.sponge

    spongealphaz = spongealphaz_dim * tref

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            spongealphaz = spongealphaz / 3.0
            spongealphax = spongealphaz
            spongealphay = spongealphaz
        elseif sizex > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = spongealphaz
            spongealphay = 0.0
        elseif sizey > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = 0.0
            spongealphay = spongealphaz
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphaunifiedsponge .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphaunifiedsponge[i, j, k] =
                    alphaunifiedsponge[i, j, k] +
                    spongealphaz * ((height - zsponge) / dzsponge)^spongeorder
            end
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= xsponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax *
                        ((xsponge0 - x[io + i]) / dxsponge)^spongeorder
                elseif x[io + i] >= xsponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax *
                        ((x[io + i] - xsponge1) / dxsponge)^spongeorder
                end
            end
            if sizey > 1
                if y[jo + j] <= ysponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay *
                        ((ysponge0 - y[jo + j]) / dysponge)^spongeorder
                elseif y[jo + j] >= ysponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay *
                        ((y[jo + j] - ysponge1) / dysponge)^spongeorder
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphaunifiedsponge, namelists, domain)
        set_meridional_boundaries_of_field!(
            alphaunifiedsponge,
            namelists,
            domain,
        )
    end

    ko == 0 &&
        @views alphaunifiedsponge[:, :, k0 - 1] .= alphaunifiedsponge[:, :, k0]
    ko + nzz == sizezz &&
        @views alphaunifiedsponge[:, :, k1 + 1] .= alphaunifiedsponge[:, :, k1]

    return
end

function compute_sponge!(
    state::State,
    dt::AbstractFloat,
    spongetype::SinusoidalSponge,
)
    (; namelists, domain) = state
    (; sizex, sizey, sizez) = namelists.domain
    (; sizezz, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) = domain
    (; x, y, ztfc) = state.grid
    (; tref) = state.constants
    (; lateralsponge, spongealphaz_dim) = namelists.sponge
    (;
        alphaunifiedsponge,
        dxsponge,
        dysponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        zsponge,
    ) = state.sponge

    spongealphaz = spongealphaz_dim * tref

    if lateralsponge
        ix0 = i0
        ix1 = i1
        jy0 = j0
        jy1 = j1

        if sizex > 1 && sizey > 1
            spongealphaz = spongealphaz / 3.0
            spongealphax = spongealphaz
            spongealphay = spongealphaz
        elseif sizex > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = spongealphaz
            spongealphay = 0.0
        elseif sizey > 1
            spongealphaz = spongealphaz / 2.0
            spongealphax = 0.0
            spongealphay = spongealphaz
        end
    else
        ix0 = 1
        ix1 = nxx
        jy0 = 1
        jy1 = nyy
    end

    alphaunifiedsponge .= 0.0

    kz0 = ko == 0 ? k0 : k0 - 1
    kz1 = ko + nzz == sizezz ? k1 : k1 + 1

    for k in kz0:kz1, j in jy0:jy1, i in ix0:ix1
        height = ztfc[i, j, k]

        if sizez > 1
            if height >= zsponge
                alphaunifiedsponge[i, j, k] =
                    alphaunifiedsponge[i, j, k] +
                    spongealphaz *
                    sin(0.5 * pi * (height - zsponge) / dzsponge)^2.0
            end
        end
        if lateralsponge
            if sizex > 1
                if x[io + i] <= xsponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax *
                        sin(0.5 * pi * (xsponge0 - x[io + i]) / dxsponge)^2.0
                elseif x[io + i] >= xsponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphax *
                        sin(0.5 * pi * (x[io + i] - xsponge1) / dxsponge)^2.0
                end
            end
            if sizey > 1
                if y[jo + j] <= ysponge0
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay *
                        sin(0.5 * pi * (ysponge0 - y[jo + j]) / dysponge)^2.0
                elseif y[jo + j] >= ysponge1
                    alphaunifiedsponge[i, j, k] =
                        alphaunifiedsponge[i, j, k] +
                        spongealphay *
                        sin(0.5 * pi * (y[jo + j] - ysponge1) / dysponge)^2.0
                end
            end
        end
    end

    if lateralsponge
        set_zonal_boundaries_of_field!(alphaunifiedsponge, namelists, domain)
        set_meridional_boundaries_of_field!(
            alphaunifiedsponge,
            namelists,
            domain,
        )
    end

    ko == 0 &&
        @views alphaunifiedsponge[:, :, k0 - 1] .= alphaunifiedsponge[:, :, k0]
    ko + nzz == sizezz &&
        @views alphaunifiedsponge[:, :, k1 + 1] .= alphaunifiedsponge[:, :, k1]

    return
end

"""
```julia
compute_sponge!(state::State, dt::AbstractFloat)
```

Compute sponge layer damping coefficients for wave absorption.

Dispatches to unified or legacy sponge computation based on configuration.
Sponge layers prevent spurious wave reflections from domain boundaries through
gradual damping of wave amplitudes.

# Arguments

  - `state::State`: Complete simulation state
  - `dt::AbstractFloat`: Time step size for legacy damping coefficient scaling

# Sponge Layer Design

**Purpose**: Absorb outgoing waves near domain boundaries without artificial reflections

**Implementation**:

  - **Unified sponge**: Single coefficient array for all variables
  - **Legacy sponge**: Separate coefficients for different grid staggerings
"""
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

"""
```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::ExponentialSponge)
```

Compute exponential sponge layer with distance-based damping.

Implements damping coefficient: `α(r) = α₀ * exp((r - r_boundary) / δ)`
where δ is the decay length scale.

# Features

  - **Smooth onset**: Exponential transition prevents discontinuities
  - **Physical realism**: Mimics atmospheric dissipation processes
  - **Configurable decay**: Different length scales for each direction
"""
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

"""
```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::COSMOSponge)
```

Compute COSMO-style sponge layer using cosine profile.

Implements damping coefficient: `α(r) = α₀/2 * (1 - cos(π * (r - r_start) / δ))`.
"""
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

"""
```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::PolynomialSponge)
```

Compute polynomial sponge layer with power-law damping profile.

Implements damping coefficient: `α(r) = α₀ * ((r - r_start) / δ)^n`
where n is the polynomial order.

# Features

  - **Flexible profiles**: Adjustable polynomial order for different transition shapes
  - **Sharp transitions**: Higher orders create more localized damping regions
  - **Computational efficiency**: Simple power function evaluation
"""
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

"""
```julia
compute_sponge!(state::State, dt::AbstractFloat, spongetype::SinusoidalSponge)
```

Compute sinusoidal sponge layer using squared sine profile.

Implements damping coefficient: `α(r) = α₀ * sin²(π/2 * (r - r_start) / δ)`
providing smooth quadratic-like transitions.

# Features

  - **Smooth transitions**: sin² profile ensures smooth derivatives
  - **Moderate gradients**: Balanced between sharp and gradual transitions
"""
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

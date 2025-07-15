"""
```julia
Predictands{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Storage for prognostic variables (predictands) on a 3D grid.

# Fields

  - `rho::A`: Density field (nxx × nyy × nzz)
  - `rhop::A`: Density perturbation field (nxx × nyy × nzz)
  - `u::A`: x-velocity field (nxx × nyy × nzz)
  - `v::A`: y-velocity field (nxx × nyy × nzz)
  - `w::A`: z-velocity field (nxx × nyy × nzz)
  - `pip::A`: Pressure perturbation field (nxx × nyy × nzz)
  - `p::B`: Pressure field (model-dependent dimensions)

For incompressible models, `p` is empty (0×0×0).
"""
struct Predictands{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    rho::A
    rhop::A
    u::A
    v::A
    w::A
    pip::A
    p::B
end

"""
```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)
```

Construct `Predictands` from configuration namelists, constants, domain, and atmosphere.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)
    (; model, testcase) = namelists.setting
    return Predictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        model,
        testcase,
    )
end

"""
```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    model::AbstractModel,
    testcase::AbstractTestCase,
)
```

Construct `Predictands` for incompressible models. Initializes velocity fields with background flow and sets empty pressure field.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::AbstractModel,
    testcase::AbstractTestCase,
)
    (; backgroundflow_dim) = namelists.atmosphere
    (; uref) = constants
    (; nxx, nyy, nzz) = domain

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    p = zeros(0, 0, 0)

    # Set the initial winds.
    u .= backgroundflow_dim[1] / uref
    v .= backgroundflow_dim[2] / uref
    w .= backgroundflow_dim[3] / uref

    # Return a Predictands instance.
    return Predictands(rho, rhop, u, v, w, pip, p)
end

"""
```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    model::Compressible,
    testcase::AbstractTestCase,
)
```

Construct `Predictands` for compressible models. Initializes velocity fields with background flow and pressure field with stratified values.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::Compressible,
    testcase::AbstractTestCase,
)
    (; backgroundflow_dim) = namelists.atmosphere
    (; uref) = constants
    (; nxx, nyy, nzz) = domain
    (; pstrattfc) = atmosphere

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip, p) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Set the initial winds.
    u .= backgroundflow_dim[1] / uref
    v .= backgroundflow_dim[2] / uref
    w .= backgroundflow_dim[3] / uref

    # Set the initial mass-weighted potential temperature.
    p .= pstrattfc

    # Return a Predictands instance.
    return Predictands(rho, rhop, u, v, w, pip, p)
end

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    model::PseudoIncompressible,
    testcase::WavePacket,
)
    (; backgroundflow_dim, theta0_dim) = namelists.atmosphere
    (; nbx, nby, nbz) = namelists.domain
    (;
        lambdax_dim,
        lambday_dim,
        lambdaz_dim,
        x0_dim,
        y0_dim,
        z0_dim,
        sigmax_dim,
        sigmay_dim,
        sigmaz_dim,
        a0,
        branch,
        wavepacketdim,
    ) = namelists.wavepacket
    (; uref, lref, tref, kappa, ma, thetaref, fr2) = constants
    (; nxx, nyy, nzz, k0, k1, j0, j1, i0, i1, io, jo) = domain
    (; bvsstrattfc, fc, rhostrattfc, pstrattfc) = atmosphere
    (; x, y, ztfc, jac, met) = grid

    # Initialize the predictands.
    (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    p = zeros(0, 0, 0)

    # Set the initial winds.
    u .= backgroundflow_dim[1] / uref
    v .= backgroundflow_dim[2] / uref
    w .= backgroundflow_dim[3] / uref

    if lambdax_dim == 0.0
        kk = 0.0
    else
        kk = 2.0 * pi / (lambdax_dim / lref)
    end

    if lambday_dim == 0.0
        ll = 0.0
    else
        ll = 2.0 * pi / (lambday_dim / lref)
    end

    mm = 2.0 * pi / (lambdaz_dim / lref)

    x0 = x0_dim / lref
    y0 = y0_dim / lref
    z0 = z0_dim / lref

    sigmax = sigmax_dim / lref
    sigmay = sigmay_dim / lref
    sigmaz = sigmaz_dim / lref

    kh = sqrt(kk^2.0 + ll^2.0)

    for kz in (k0 - nbz):(k1 + nbz),
        jy in (j0 - nby):(j1 + nby),
        ix in (i0 - nbx):(i1 + nbx)

        n2 = bvsstrattfc[ix, jy, kz]
        f = fc[jy]
        f2 = f^2.0
        omega = branch * sqrt((n2 * kh^2.0 + f2 * mm^2.0) / (kh^2.0 + mm^2.0))

        omega2 = omega^2.0

        if wavepacketdim == 1
            deltax = 0.0
            deltay = 0.0
        elseif wavepacketdim == 2
            deltax = x[io + ix] - x0
            deltay = 0.0
        elseif wavepacketdim == 3
            deltax = (x[io + ix] - x0)
            deltay = (y[jo + jy] - y0)
        end

        deltaz = ztfc[ix, jy, kz] - z0

        # Gaussian in z, Cosine in x and y.

        if sigmax == 0.0
            envel = 1.0
        elseif abs(deltax) < sigmax
            envel = 0.5 * (1.0 + cos(deltax * pi / sigmax))
        else
            envel = 0.0
        end

        if sigmay == 0.0
            envel = 1.0 * envel
        elseif abs(deltay) < sigmay
            envel = envel * 0.5 * (1.0 + cos(deltay * pi / sigmay))
        else
            envel = 0.0
        end

        envel = envel * exp(-deltaz^2.0 / (2.0 * sigmaz^2.0))

        theta0 = theta0_dim / thetaref

        bhat = a0 * n2 / mm * envel
        uhat =
            1im / (mm * n2) * (omega2 - n2) / (omega2 - f2) *
            (kk * omega + 1im * ll * f) *
            bhat
        vhat =
            1im / (mm * n2) * (omega2 - n2) / (omega2 - f2) *
            (ll * omega - 1im * kk * f) *
            bhat
        what = 1im * omega / n2 * bhat
        piphat = 1im * kappa * ma^2.0 * (omega2 - n2) / n2 / mm / theta0 * bhat

        phi = kk * x[io + ix] + ll * y[jo + jy] + mm * ztfc[ix, jy, kz]

        bprime = real(bhat * exp(1im * phi))
        uprime = real(uhat * exp(1im * phi))
        vprime = real(vhat * exp(1im * phi))
        wprime = real(what * exp(1im * phi))
        pipprime = real(piphat * exp(1im * phi))

        rhoprime = 1.0 / (1.0 + fr2 * bprime) * rhostrattfc[ix, jy, kz]
        rhoprime = rhoprime - rhostrattfc[ix, jy, kz]

        u[ix, jy, kz] = u[ix, jy, kz] + uprime
        v[ix, jy, kz] = v[ix, jy, kz] + vprime
        w[ix, jy, kz] = w[ix, jy, kz] + wprime
        rho[ix, jy, kz] = rhoprime
        pip[ix, jy, kz] = pipprime

        w[ix, jy, kz] =
            w[ix, jy, kz] / jac[ix, jy, kz] +
            met[ix, jy, kz, 1, 3] * u[ix, jy, kz] +
            met[ix, jy, kz, 2, 3] * v[ix, jy, kz]
    end

    for ix in (i0 - 1):(i1 + 1)
        u[ix, :, :] .= 0.5 .* (u[ix, :, :] .+ u[ix + 1, :, :])
    end

    for jy in (j0 - 1):(j1 + 1)
        v[:, jy, :] .= 0.5 .* (v[:, jy, :] .+ v[:, jy + 1, :])
    end

    for kz in (k0 - 1):(k1 + 1),
        jy in (j0 - nby):(j1 + nby),
        ix in (i0 - nbx):(i1 + nbx)

        w[ix, jy, kz] =
            (
                jac[ix, jy, kz + 1] * w[ix, jy, kz] +
                jac[ix, jy, kz] * w[ix, jy, kz + 1]
            ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
    end

    w[:, :, k0] .= 0.0
    w[:, :, k1] .= 0.0

    rhop .= rho

    # Return a Predictands instance.
    return Predictands(rho, rhop, u, v, w, pip, p)
end

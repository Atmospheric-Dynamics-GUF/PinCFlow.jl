"""
```julia
Predictands{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
```

Arrays for prognostic variables.

```julia
Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)::Predictands
```

Construct a `Predictands` instance.

The predictands are initialized with the corresponding functions in `namelists.atmosphere`. The mass-weighted potential temperature `p` is constructed depending on the dynamic equations (see `set_p`).

# Fields

  - `rho::A`: Density.

  - `rhop::A`: Density-fluctuations.

  - `u::A`: Zonal wind.

  - `v::A`: Meridional wind.

  - `w::A`: Transformed vertical wind.

  - `pip::A`: Exner-pressure fluctuations.

  - `p::B`: Mass-weighted potential temperature.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields that describe the grid.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.VariableTypes.set_p`](@ref)
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

function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)::Predictands
    (;
        initial_rhop,
        initial_thetap,
        initial_u,
        initial_v,
        initial_w,
        initial_pip,
    ) = namelists.atmosphere
    (; model) = namelists.atmosphere
    (; lref, rhoref, thetaref, uref) = constants
    (; i0, i1, j0, j1, k0, k1, nxx, nyy, nzz) = domain
    (; x, y, zc, met, jac) = grid
    (; rhobar, thetabar) = atmosphere

    (rho, rhop, thetap, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:7)

    @ivy for k in 1:nzz, j in j0:j1, i in i0:i1
        xdim = x[i] * lref
        ydim = y[j] * lref
        zcdim = zc[i, j, k] * lref

        rhop[i, j, k] = initial_rhop(xdim, ydim, zcdim) / rhoref
        thetap[i, j, k] = initial_thetap(xdim, ydim, zcdim) / thetaref
        u[i, j, k] = initial_u(xdim, ydim, zcdim) / uref
        v[i, j, k] = initial_v(xdim, ydim, zcdim) / uref
        w[i, j, k] = initial_w(xdim, ydim, zcdim) / uref
        pip[i, j, k] = initial_pip(xdim, ydim, zcdim)
    end

    for f! in
        (set_zonal_boundaries_of_field!, set_meridional_boundaries_of_field!)
        f!(rhop, namelists, domain)
        f!(thetap, namelists, domain)
        f!(u, namelists, domain)
        f!(v, namelists, domain)
        f!(w, namelists, domain)
        f!(pip, namelists, domain)
    end

    rho .= rhop

    @ivy w .= met[:, :, :, 1, 3] .* u .+ met[:, :, :, 2, 3] .* v .+ w ./ jac

    @ivy for i in i0:i1
        u[i, :, :] .= (u[i, :, :] .+ u[i + 1, :, :]) ./ 2
    end
    set_zonal_boundaries_of_field!(u, namelists, domain)

    @ivy for j in j0:j1
        v[:, j, :] .= (v[:, j, :] .+ v[:, j + 1, :]) ./ 2
    end
    set_meridional_boundaries_of_field!(v, namelists, domain)

    @ivy for k in k0:k1
        w[:, :, k] .=
            (
                jac[:, :, k + 1] .* w[:, :, k] .+
                jac[:, :, k] .* w[:, :, k + 1]
            ) ./ (jac[:, :, k] .+ jac[:, :, k + 1])
    end
    set_vertical_boundaries_of_field!(w, namelists, domain, -; staggered = true)

    p = set_p(model, rhobar, thetabar, rhop, thetap)

    return Predictands(rho, rhop, u, v, w, pip, p)
end

"""
```julia
Constants{A <: AbstractFloat}
```

Physical constants and reference quantities.

This struct contains all dimensional and non-dimensional parameters used throughout the model, including thermodynamic constants, reference scales, and derived flow parameters for proper non-dimensionalization.

# Fields

## Physical Constants

  - `gamma::A`: Ratio of specific heats ``c_p/c_v = 1.4``
  - `gammainv::A`: Inverse ratio of specific heats ``1/\\gamma``
  - `kappa::A`: Poisson constant ``(\\gamma-1)/\\gamma = R/c_p ≈ 0.286``
  - `kappainv::A`: Inverse Poisson constant ``\\gamma/(\\gamma-1) = c_p/R``
  - `rsp::A`: Specific gas constant ``R`` [J/(kg·K)] = 287.0
  - `g::A`: Gravitational acceleration [m/s²] = 9.81

## Reference Quantities

  - `rhoref::A`: Reference density [kg/m³] ≡ 1.184 (sea level)
  - `pref::A`: Reference pressure [Pa] ≡ 101325.0 (sea level)
  - `aref::A`: Reference sound speed [m/s] ``= \\sqrt{p_\\mathrm{ref}/\\rho_\\mathrm{ref}}``
  - `uref::A`: Reference velocity [m/s] ``= a_\\mathrm{ref}`` (Mach 1 scaling)
  - `lref::A`: Reference length [m] ``= p_\\mathrm{ref}/(\\rho_\\mathrm{ref} \\cdot g)`` (pressure scale height)
  - `tref::A`: Reference time [s] ``= l_\\mathrm{ref}/a_\\mathrm{ref}`` (acoustic time scale)
  - `thetaref::A`: Reference temperature [K] ``= a_\\mathrm{ref}^2/r_\\mathrm{sp}``
  - `fref::A`: Reference body force [N/m³] ``= \\rho_\\mathrm{ref} \\cdot u_\\mathrm{ref}^2/l_\\mathrm{ref}``

## Non-dimensional Parameters

  - `g_ndim::A`: Non-dimensional gravity ``= g/(u_\\mathrm{ref}^2/l_\\mathrm{ref}) ≡ \\mathrm{Fr}^{-2}``
  - `re::A`: Reynolds number ``= \\rho_\\mathrm{ref} \\cdot u_\\mathrm{ref} \\cdot l_\\mathrm{ref}/\\mu``
  - `ma::A`: Mach number ``= u_\\mathrm{ref}/a_\\mathrm{ref} ≡ 1.0`` (by construction)
  - `mainv2::A`: Inverse Mach number squared ``= (a_\\mathrm{ref}/u_\\mathrm{ref})^2``
  - `ma2::A`: Mach number squared ``= (u_\\mathrm{ref}/a_\\mathrm{ref})^2``
  - `fr::A`: Froude number ``= u_\\mathrm{ref}/\\sqrt{g \\cdot l_\\mathrm{ref}} ≡ 1.0`` (by construction)
  - `frinv2::A`: Inverse Froude number squared ``= g \\cdot l_\\mathrm{ref}/u_\\mathrm{ref}^2``
  - `fr2::A`: Froude number squared ``= u_\\mathrm{ref}^2/(g \\cdot l_\\mathrm{ref})``
  - `sig::A`: Stratification parameter ``= \\mathrm{Ma}^2/\\mathrm{Fr}^2 = \\rho_\\mathrm{ref} \\cdot g \\cdot l_\\mathrm{ref}/p_\\mathrm{ref}``
"""
struct Constants{A <: AbstractFloat}

    # Natural constants.
    gamma::A
    gammainv::A
    kappa::A
    kappainv::A
    rsp::A
    g::A

    # Reference quantities.
    rhoref::A
    pref::A
    aref::A
    uref::A
    lref::A
    tref::A
    thetaref::A
    fref::A

    # Non-dimensionalized gravitational acceleration.
    g_ndim::A

    # Flow parameters.
    re::A
    ma::A
    mainv2::A
    ma2::A
    fr::A
    frinv2::A
    fr2::A
    sig::A
end

"""
```julia
Constants(namelists::Namelists)
```

Creates a `Constants` instance from simulation parameters.

# Arguments

  - `namelists`: Configuration containing:

      + `atmosphere.specifyreynolds`: Whether to use prescribed Reynolds number
      + `atmosphere.reinv`: Inverse Reynolds number (if prescribed)
      + `atmosphere.mu_viscous_dim`: Dynamic viscosity (if computed)

# Non-dimensionalization scheme

The model uses a non-dimensionalization based on:

  - **Length scale**: Pressure scale height ``l_\\mathrm{ref} = p_\\mathrm{ref}/(\\rho_\\mathrm{ref}g)``
  - **Velocity scale**: Sound speed ``u_\\mathrm{ref} = \\sqrt{p_\\mathrm{ref}/\\rho_\\mathrm{ref}}``
  - **Time scale**: Acoustic time ``t_\\mathrm{ref} = l_\\mathrm{ref}/u_\\mathrm{ref}``
  - **Pressure scale**: Reference pressure ``p_\\mathrm{ref}``
  - **Density scale**: Reference density ``\\rho_\\mathrm{ref}``
  - **Temperature scale**: Acoustic temperature ``\\theta_\\mathrm{ref} = u_\\mathrm{ref}^2/r_\\mathrm{sp}``
"""
function Constants(namelists::Namelists)

    # Get parameters.
    (; specifyreynolds, reinv, mu_viscous_dim) = namelists.atmosphere

    # Set natural constants.
    gamma = 1.4
    gammainv = 1.0 / gamma
    kappa = (gamma - 1.0) / gamma
    kappainv = 1.0 / kappa
    rsp = 287.0
    g = 9.81

    # Set reference quantites.
    rhoref = 1.184 # in kg/m^3
    pref = 101325.0 # in Pa = kg/m/s^2
    aref = sqrt(pref / rhoref) # in m/s
    uref = aref # in m/s
    lref = pref / rhoref / g # in m
    tref = lref / aref # in s
    thetaref = aref^2 / rsp # in K
    fref = rhoref * uref^2 / lref # in N/m^3

    # Compute non-dimensionalized gravitational acceleration.
    g_ndim = g / (uref^2 / lref)

    # Set the Reynolds number.
    if specifyreynolds
        if reinv < eps()
            re = 1 / eps()
        else
            re = 1.0 / reinv
        end
    else
        if mu_viscous_dim / uref / lref < eps()
            re = 1 / eps()
        else
            re = uref * lref / mu_viscous_dim
        end
    end

    # Set other flow parameters.
    ma = uref / aref # Ma = 1
    mainv2 = 1.0 / ma^2
    ma2 = ma^2
    fr = uref / sqrt(g * lref) # Fr = 1
    frinv2 = 1.0 / fr^2
    fr2 = fr^2
    sig = ma^2 / fr^2

    # Return a Constants instance.
    return Constants(
        gamma,
        gammainv,
        kappa,
        kappainv,
        rsp,
        g,
        rhoref,
        pref,
        aref,
        uref,
        lref,
        tref,
        thetaref,
        fref,
        g_ndim,
        re,
        ma,
        mainv2,
        ma2,
        fr,
        frinv2,
        fr2,
        sig,
    )
end

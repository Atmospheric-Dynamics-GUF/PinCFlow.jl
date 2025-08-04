"""
```julia
Constants{A <: AbstractFloat}
```

Composite type for natural constants, reference quantities and non-dimensional parameters.

# Fields

Natural constants:

  - `gamma::A`: Ratio of specific heats ``\\gamma = c_p / c_V = 1.4``.
  - `gammainv::A`: Inverse ratio of specific heats ``1 / \\gamma``.
  - `kappa::A`: Ratio between specific gas constant and specific heat capacity at constant pressure ``\\kappa = \\left(\\gamma - 1\\right) / \\gamma = R / c_p = 2 / 7``.
  - `kappainv::A`: Ratio between specific heat capacity at constant pressure and specific gas constant ``1 / \\kappa``.
  - `rsp::A`: Specific gas constant ``R = 287 \\, \\mathrm{J \\, kg^{- 1} \\, K^{- 1}}``.
  - `g::A`: Gravitational acceleration ``g = 9.81 \\, \\mathrm{m \\, s^{- 2}}``.

Reference quantities:

  - `rhoref::A`: Reference density ``\\rho_\\mathrm{ref} = 1.184 \\, \\mathrm{kg \\, m^{- 3}}``.
  - `pref::A`: Reference pressure ``p_\\mathrm{ref} = 101325 \\, \\mathrm{Pa}``.
  - `aref::A`: Reference sound speed ``c_\\mathrm{ref} = \\sqrt{p_\\mathrm{ref} / \\rho_\\mathrm{ref}}``.
  - `uref::A`: Reference wind ``u_\\mathrm{ref} = a_\\mathrm{ref}``.
  - `lref::A`: Reference length ``L_\\mathrm{ref} = p_\\mathrm{ref} /\\left(g \\rho_\\mathrm{ref}\\right)``.
  - `tref::A`: Reference time ``t_\\mathrm{ref} = L_\\mathrm{ref} / a_\\mathrm{ref}``.
  - `thetaref::A`: Reference potential temperature ``\\theta_\\mathrm{ref} = a_\\mathrm{ref}^2 / R``.
  - `fref::A`: Reference body force ``F_\\mathrm{ref} = \\rho_\\mathrm{ref} u_\\mathrm{ref}^2 / L_\\mathrm{ref}``.

Non-dimensional parameters

  - `g_ndim::A`: Non-dimensional gravitational acceleration ``\\widehat{g} = g L_\\mathrm{ref} / u_\\mathrm{ref}^2``.
  - `re::A`: Reynolds number ``\\mathrm{Re} = L_\\mathrm{ref} u_\\mathrm{ref} / \\mu`` (with ``\\mu`` being the kinematic viscosity at the surface).
  - `ma::A`: Mach number ``\\mathrm{Ma} = u_\\mathrm{ref} / a_\\mathrm{ref}``.
  - `mainv2::A`: Inverse Mach number squared ``\\mathrm{Ma}^{- 2}``.
  - `ma2::A`: Mach number squared ``\\mathrm{Ma}^2``.
  - `fr::A`: Froude number ``\\mathrm{Fr} = u_\\mathrm{ref} / \\sqrt{g L_\\mathrm{ref}}``.
  - `frinv2::A`: Inverse Froude number squared ``\\mathrm{Fr}^{- 2}``.
  - `fr2::A`: Froude number squared ``\\mathrm{Fr}^{2}``.
  - `sig::A`: Ratio between squared Mach number and squared Froude number ``\\sigma = \\mathrm{Ma}^2 / \\mathrm{Fr}^2``.
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

Creates a `Constants` instance.

The Reynolds number is the only constant that depends on the model parameters in `namelists`. If `namelists.atmosphere.specifyreynolds` is `false`, the Reynolds number is ``\\mathrm{Re} = L_\\mathrm{ref} u_\\mathrm{ref} / \\mu``, with ``\\mu`` being the kinematic viscosity at the surface, given by `namelists.atmosphere.mu_viscous_dim`. Otherwise, it is set to the inverse of `namelists.atmosphere.reinv`.

# Arguments

  - `namelists`: Namelists with all model parameters.

# Returns

  - `::Constants`: `Constants` instance.
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

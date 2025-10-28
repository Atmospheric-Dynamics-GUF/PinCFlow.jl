"""
```julia
Constants{
    A <: AbstractFloat,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: AbstractFloat,
    M <: AbstractFloat,
    N <: AbstractFloat,
    O <: AbstractFloat,
    P <: AbstractFloat,
    Q <: AbstractFloat,
    R <: AbstractFloat,
    S <: AbstractFloat,
    T <: AbstractFloat,
    U <: AbstractFloat,
    V <: AbstractFloat,
    W <: AbstractFloat,
}
```

Composite type for natural constants, reference quantities and non-dimensional parameters.

```julia
Constants(namelists::Namelists)::Constants
```

Create a `Constants` instance.

The Reynolds number is the only constant that depends on the model parameters in `namelists`. If `namelists.atmosphere.specify_reynolds_number` is `false`, the Reynolds number is ``\\mathrm{Re} = L_\\mathrm{ref} u_\\mathrm{ref} / \\mu``, with ``\\mu`` being the kinematic viscosity at the surface, given by `namelists.atmosphere.kinematic_viscosity`. Otherwise, it is set to the inverse of `namelists.atmosphere.inverse_reynolds_number`.

# Fields

Natural constants:

  - `gamma::A`: Ratio of specific heats ``\\gamma = c_p / c_V = 1.4``.

  - `gammainv::B`: Inverse ratio of specific heats ``1 / \\gamma``.

  - `kappa::C`: Ratio between specific gas constant and specific heat capacity at constant pressure ``\\kappa = \\left(\\gamma - 1\\right) / \\gamma = R / c_p = 2 / 7``.

  - `kappainv::D`: Ratio between specific heat capacity at constant pressure and specific gas constant ``1 / \\kappa``.

  - `rsp::E`: Specific gas constant ``R = 287 \\, \\mathrm{J \\, kg^{- 1} \\, K^{- 1}}``.

  - `g::F`: Gravitational acceleration ``g = 9.81 \\, \\mathrm{m \\, s^{- 2}}``.

Reference quantities:

  - `rhoref::G`: Reference density ``\\rho_\\mathrm{ref} = 1.184 \\, \\mathrm{kg \\, m^{- 3}}``.

  - `pref::H`: Reference pressure ``p_\\mathrm{ref} = 101325 \\, \\mathrm{Pa}``.

  - `aref::I`: Reference sound speed ``c_\\mathrm{ref} = \\sqrt{p_\\mathrm{ref} / \\rho_\\mathrm{ref}}``.

  - `uref::J`: Reference wind ``u_\\mathrm{ref} = a_\\mathrm{ref}``.

  - `lref::K`: Reference length ``L_\\mathrm{ref} = p_\\mathrm{ref} /\\left(g \\rho_\\mathrm{ref}\\right)``.

  - `tref::L`: Reference time ``t_\\mathrm{ref} = L_\\mathrm{ref} / a_\\mathrm{ref}``.

  - `thetaref::M`: Reference potential temperature ``\\theta_\\mathrm{ref} = a_\\mathrm{ref}^2 / R``.

  - `fref::N`: Reference body force ``F_\\mathrm{ref} = \\rho_\\mathrm{ref} u_\\mathrm{ref}^2 / L_\\mathrm{ref}``.

Non-dimensional parameters

  - `g_ndim::O`: Non-dimensional gravitational acceleration ``\\widehat{g} = g L_\\mathrm{ref} / u_\\mathrm{ref}^2``.

  - `re::P`: Reynolds number ``\\mathrm{Re} = L_\\mathrm{ref} u_\\mathrm{ref} / \\mu`` (with ``\\mu`` being the kinematic viscosity at the surface).

  - `ma::Q`: Mach number ``\\mathrm{Ma} = u_\\mathrm{ref} / a_\\mathrm{ref}``.

  - `mainv2::R`: Inverse Mach number squared ``\\mathrm{Ma}^{- 2}``.

  - `ma2::S`: Mach number squared ``\\mathrm{Ma}^2``.

  - `fr::T`: Froude number ``\\mathrm{Fr} = u_\\mathrm{ref} / \\sqrt{g L_\\mathrm{ref}}``.

  - `frinv2::U`: Inverse Froude number squared ``\\mathrm{Fr}^{- 2}``.

  - `fr2::V`: Froude number squared ``\\mathrm{Fr}^{2}``.

  - `sig::W`: Ratio between squared Mach number and squared Froude number ``\\sigma = \\mathrm{Ma}^2 / \\mathrm{Fr}^2``.

# Arguments

  - `namelists`: Namelists with all model parameters.
"""
struct Constants{
    A <: AbstractFloat,
    B <: AbstractFloat,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: AbstractFloat,
    M <: AbstractFloat,
    N <: AbstractFloat,
    O <: AbstractFloat,
    P <: AbstractFloat,
    Q <: AbstractFloat,
    R <: AbstractFloat,
    S <: AbstractFloat,
    T <: AbstractFloat,
    U <: AbstractFloat,
    V <: AbstractFloat,
    W <: AbstractFloat,
}

    # Natural constants.
    gamma::A
    gammainv::B
    kappa::C
    kappainv::D
    rsp::E
    g::F

    # Reference quantities.
    rhoref::G
    pref::H
    aref::I
    uref::J
    lref::K
    tref::L
    thetaref::M
    fref::N

    # Non-dimensionalized gravitational acceleration.
    g_ndim::O

    # Flow parameters.
    re::P
    ma::Q
    mainv2::R
    ma2::S
    fr::T
    frinv2::U
    fr2::V
    sig::W
end

function Constants(namelists::Namelists)::Constants
    (; specify_reynolds_number, inverse_reynolds_number, kinematic_viscosity) =
        namelists.atmosphere

    # Set natural constants.
    gamma = 1.4
    gammainv = 1.0 / gamma
    kappa = (gamma - 1.0) / gamma
    kappainv = 1.0 / kappa
    rsp = 287.0
    g = 9.81

    # Set reference quantities.
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
    if specify_reynolds_number
        if inverse_reynolds_number < eps()
            re = 1 / eps()
        else
            re = 1.0 / inverse_reynolds_number
        end
    else
        if kinematic_viscosity / uref / lref < eps()
            re = 1 / eps()
        else
            re = uref * lref / kinematic_viscosity
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

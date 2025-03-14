
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

  # Small float.
  small::A
end

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
    if reinv < 1.0e-20
      re = 1.0e20
    else
      re = 1.0 / reinv
    end
  else
    if mu_viscous_dim / uref / lref < 1.0e-20
      re = 1.0e20
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

  # Set small float.
  # small = 1.0e-20
  small = nextfloat(0.0)

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
    small,
  )
end

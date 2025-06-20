"""
    Constants{A<:AbstractFloat}

Physical constants and reference quantities.

This struct contains all dimensional and non-dimensional parameters used throughout
the model, including thermodynamic constants, reference scales, and
derived flow parameters for proper non-dimensionalization.

# Type Parameters

  - `A<:AbstractFloat`: Floating-point precision type

# Physical Constants

  - `gamma::A`: Ratio of specific heats for dry air (cp/cv ≡ 1.4)
  - `gammainv::A`: Inverse ratio of specific heats (1/γ)
  - `kappa::A`: Poisson constant (γ-1)/γ ≡ R/cp ≈ 0.286
  - `kappainv::A`: Inverse Poisson constant γ/(γ-1) ≡ cp/R
  - `rsp::A`: Specific gas constant for dry air [J/(kg·K)] ≡ 287.0
  - `g::A`: Gravitational acceleration [m/s²] ≡ 9.81

# Reference Quantities

  - `rhoref::A`: Reference density [kg/m³] ≡ 1.184 (sea level)
  - `pref::A`: Reference pressure [Pa] ≡ 101325.0 (sea level)
  - `aref::A`: Reference sound speed [m/s] = √(pref/ρref)
  - `uref::A`: Reference velocity [m/s] = aref (Mach 1 scaling)
  - `lref::A`: Reference length [m] = pref/(ρref·g) (pressure scale height)
  - `tref::A`: Reference time [s] = lref/aref (acoustic time scale)
  - `thetaref::A`: Reference temperature [K] = aref²/rsp
  - `fref::A`: Reference body force [N/m³] = ρref·uref²/lref

# Non-dimensional Parameters

  - `g_ndim::A`: Non-dimensional gravity = g/(uref²/lref) ≡ Fr⁻²
  - `re::A`: Reynolds number = ρref·uref·lref/μ
  - `ma::A`: Mach number = uref/aref ≡ 1.0 (by construction)
  - `mainv2::A`: Inverse Mach number squared = (aref/uref)²
  - `ma2::A`: Mach number squared = (uref/aref)²
  - `fr::A`: Froude number = uref/√(g·lref) ≡ 1.0 (by construction)
  - `frinv2::A`: Inverse Froude number squared = g·lref/uref²
  - `fr2::A`: Froude number squared = uref²/(g·lref)
  - `sig::A`: Stratification parameter = Ma²/Fr² = ρref·g·lref/pref

# Constructor

    Constants(namelists::Namelists)

Creates a Constants instance from simulation parameters.

## Parameters

  - `namelists::Namelists`: Configuration containing:

      + `atmosphere.specifyreynolds`: Whether to use prescribed Reynolds number
      + `atmosphere.reinv`: Inverse Reynolds number (if prescribed)
      + `atmosphere.mu_viscous_dim`: Dynamic viscosity [Pa·s] (if computed)

## Non-dimensionalization Scheme

The model uses a consistent non-dimensionalization based on:

  - **Length scale**: Pressure scale height lref = pref/(ρref·g)
  - **Velocity scale**: Sound speed uref = √(pref/ρref)
  - **Time scale**: Acoustic time tref = lref/uref
  - **Pressure scale**: Reference pressure pref
  - **Density scale**: Reference density ρref
  - **Temperature scale**: Acoustic temperature θref = uref²/rsp

## Examples

```julia
# Create constants from namelists
constants = Constants(namelists)

# Access physical constants
γ = constants.gamma          # 1.4
R = constants.rsp           # 287.0 J/(kg·K)
g = constants.g             # 9.81 m/s²

# Reference scales
ρ₀ = constants.rhoref       # 1.184 kg/m³
p₀ = constants.pref         # 101325.0 Pa
c₀ = constants.aref         # ~293 m/s
L₀ = constants.lref         # ~8700 m

# Non-dimensional parameters
Re = constants.re           # Reynolds number
Ma = constants.ma           # Mach number (≡ 1)
Fr = constants.fr           # Froude number (≡ 1)
```
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

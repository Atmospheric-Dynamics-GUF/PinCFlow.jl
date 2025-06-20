"""
    Predictands{A, B}

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
    Predictands(namelists::Namelists, constants::Constants, domain::Domain, atmosphere::Atmosphere)

Construct `Predictands` from configuration namelists, constants, domain, and atmosphere.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)
    (; model, testcase) = namelists.setting
    return Predictands(
        namelists,
        constants,
        domain,
        atmosphere,
        model,
        testcase,
    )
end

"""
    Predictands(namelists, constants, domain, atmosphere, model::AbstractModel, testcase)

Construct `Predictands` for incompressible models. Initializes velocity fields with background flow and sets empty pressure field.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
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
    Predictands(namelists, constants, domain, atmosphere, model::Compressible, testcase)

Construct `Predictands` for compressible models. Initializes velocity fields with background flow and pressure field with stratified values.
"""
function Predictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
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

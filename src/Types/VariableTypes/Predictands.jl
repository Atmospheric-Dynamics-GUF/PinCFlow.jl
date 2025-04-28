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

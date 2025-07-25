struct SgsPredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    n::A # ice crystal number concentration n
    q::A # ice mixing ratio q
    qv::A # vapor mixing ratio qv
end

function SgsPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
    iceconstants::IceConstants
)
    (; icesetup) = namelists.ice

    return SgsPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        icesetup,
        variables,
        iceconstants
    )
end

function SgsPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::NoIce,
    variables::Variables,
    iceconstants::IceConstants 
)
    n = zeros(0, 0, 0)
    q = zeros(0, 0, 0)
    qv = zeros(0, 0, 0)

    return IcePredictands(n, q, qv)
end

function SgsPredictands(
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::AbstractIce,
    variables::Variables,
    iceconstants::IceConstants,
    icepredictands::IcePredictands

)

    (; nxx, nyy, nzz) = domain
    (; i0, i1, j0, j1, k0, k1) = domain 
    (; ztfc) = grid
    (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = atmosphere
    (; rho, rhop, u, v, w, pip, p) = variables.predictands
    (; kappainv, pref, gamma, lref) = constants
    (; press0_dim) = namelists.atmosphere
    (;epsil0hat, meanMassIce, mRef) = iceconstants 


    (; n, q, qv) = icepredictands

    sgs_n = n
    sgs_q = q
    sgs_qv = qv

    return SgsPredictands(sgs_n, sgs_q, sgs_qv)
end

function SgsPredictands(
    icepredictands::IcePredictands

)
    (; n, q, qv) = icepredictands

    sgs_n = n
    sgs_q = q
    sgs_qv = qv

    return SgsPredictands(sgs_n, sgs_q, sgs_qv)
end
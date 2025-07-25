struct SgsTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dn::A
    dq::A
    dqv::A
end

function SgsTendencies(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice
    return SgsTendencies(domain, icesetup)
end

function SgsTendencies(domain::Domain, icesetup::NoIce)
    dn = zeros(0, 0, 0)
    dq = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)

    return SgsTendencies(dn, dq, dqv)
end

function SgsTendencies(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    dn = zeros(nxx, nyy, nzz)
    dq = zeros(nxx, nyy, nzz)
    dqv = zeros(nxx, nyy, nzz)

    return SgsTendencies(dn, dq, dqv)
end

struct GW{A <: AbstractArray{<:AbstractFloat, 3}}
    wwp::A
    epp::A
    thp::A
end

function GW(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice
    return GW(domain, icesetup)
end

function GW(domain::Domain, icesetup::NoIce)
    wwp = zeros(0, 0, 0)
    epp = zeros(0, 0, 0)
    thp = zeros(0, 0, 0)

    return GW(wwp, epp, thp)
end

function GW(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    wwp = zeros(nxx, nyy, nzz)
    epp = zeros(nxx, nyy, nzz)
    thp = zeros(nxx, nyy, nzz)

    return GW(wwp, epp, thp)
end



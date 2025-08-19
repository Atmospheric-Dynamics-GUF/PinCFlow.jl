struct SgsGW{A <: AbstractArray{<:AbstractFloat, 3}}
    wwp::A
    epp::A
    thp::A
end

function SgsGW(namelists::Namelists, domain::Domain)
    (; icesetup) = namelists.ice
    return SgsGW(domain, icesetup)
end

function SgsGW(domain::Domain, icesetup::NoIce)
    wwp = zeros(0, 0, 0)
    epp = zeros(0, 0, 0)
    thp = zeros(0, 0, 0)

    return SgsGW(wwp, epp, thp)
end

function SgsGW(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    wwp = zeros(nxx, nyy, nzz)
    epp = zeros(nxx, nyy, nzz)
    thp = zeros(nxx, nyy, nzz)

    return SgsGW(wwp, epp, thp)
end



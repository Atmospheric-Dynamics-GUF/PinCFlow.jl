struct GW{A <: AbstractArray{<:AbstractFloat, 3}}
    wwp::A
    epp::A
    thp::A
end

function GW(namelists::Namelists, domain::Domain)
    (; ice_setup) = namelists.ice
    return GW(domain, ice_setup)
end

function GW(domain::Domain, ice_setup::NoIce)

    wwp = zeros(0, 0, 0)
    epp = zeros(0, 0, 0)
    thp = zeros(0, 0, 0)

    return GW(wwp, epp, thp)
end

function GW(domain::Domain, ice_setup::AbstractIce)
    
    (; nxx, nyy, nzz) = domain

    wwp = zeros(nxx, nyy, nzz)
    epp = zeros(nxx, nyy, nzz)
    thp = zeros(nxx, nyy, nzz)

    return GW(wwp, epp, thp)
end



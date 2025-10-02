struct SgsGW{A <: AbstractArray{<:AbstractFloat, 3}}
    wwp::A
    epp::A
    thp::A
end

function SgsGW(namelists::Namelists, 
    domain :: Domain,
    subgrid::SubGrid
    )

    (; icesetup, cloudcover) = namelists.ice

    return SgsGW(namelists, domain, subgrid, icesetup, cloudcover)

end

function SgsGW(namelists :: Namelists, 
    domain :: Domain,
    subgrid ::SubGrid,
    icesetup :: NoIce,
    cloudcover:: AbstractCloudCover
    )

    wwp = zeros(0, 0, 0)
    epp = zeros(0, 0, 0)
    thp = zeros(0, 0, 0)

    return SgsGW(wwp, epp, thp)
end

function SgsGW( namelists::Namelists, 
    domain :: Domain,
    subgrid::SubGrid,
    icesetup :: IceOn,
    cloudcover :: CloudCoverOff
    )

    (;  nxx, nyy, nzz) = domain

    wwp = zeros(nxx, nyy, nzz)
    epp = zeros(nxx, nyy, nzz)
    thp = zeros(nxx, nyy, nzz)

    return SgsGW(wwp, epp, thp)
end

function SgsGW( namelists::Namelists, 
    domain :: Domain,
    subgrid::SubGrid,
    icesetup :: AbstractIce,
    cloudcover :: CloudCoverOn
    )

    (;  nxnscxx, nynscyy, nznsczz) = subgrid
    
    wwp = zeros(nxnscxx, nynscyy, nznsczz)
    epp = zeros(nxnscxx, nynscyy, nznsczz)
    thp = zeros(nxnscxx, nynscyy, nznsczz)

    return SgsGW(wwp, epp, thp)
end



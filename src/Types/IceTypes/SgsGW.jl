struct SgsGW{A <: AbstractArray{<:AbstractFloat, 3}}
    wwp::A
    epp::A
    thp::A
end

function SgsGW(namelists::Namelists, 
    domain::Domain, 
    subgrid::SubGrid
    )

    (; icesetup) = namelists.ice

    return SgsGW(namelists, domain, subgrid, icesetup)

end

function SgsGW(namelists :: Namelists, 
    domain::Domain, 
    subgrid::SubGrid,
    icesetup ::NoIce
    )

    wwp = zeros(0, 0, 0)
    epp = zeros(0, 0, 0)
    thp = zeros(0, 0, 0)

    return SgsGW(wwp, epp, thp)
end

function SgsGW( namelists::Namelists, 
    domain::Domain, 
    subgrid::SubGrid,
    icesetup::AbstractIce
    )

    (; compute_cloudcover) = namelists.ice
    (;  nxnscxx, nynscyy, nznsczz) = subgrid
    
    if compute_cloudcover == 2

        wwp = zeros(nxnscxx, nynscyy, nznsczz)
        epp = zeros(nxnscxx, nynscyy, nznsczz)
        thp = zeros(nxnscxx, nynscyy, nznsczz)

    else 

        wwp = zeros(0, 0, 0)
        epp = zeros(0, 0, 0)
        thp = zeros(0, 0, 0)
    end

    return SgsGW(wwp, epp, thp)
end



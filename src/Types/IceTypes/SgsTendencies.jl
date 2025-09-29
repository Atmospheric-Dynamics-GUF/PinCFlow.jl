struct SgsTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dn::A
    dq::A
    dqv::A
end

function SgsTendencies(namelists::Namelists, subgrid::SubGrid)
    (; cloudcover) = namelists.ice
    return SgsTendencies(namelists, subgrid, cloudcover)
end

function SgsTendencies(namelists::Namelists, subgrid::SubGrid, cloudcover::CloudCoverOff)

    dn = zeros(0, 0, 0)
    dq = zeros(0, 0, 0)
    dqv = zeros(0, 0, 0)

    return SgsTendencies(dn, dq, dqv)
end

function SgsTendencies(namelists::Namelists, subgrid::SubGrid, cloudcover::CloudCoverOn)

    (; nxnscxx, nynscyy, nznsczz) = subgrid

    dn = zeros(nxnscxx, nynscyy, nznsczz)
    dq = zeros(nxnscxx, nynscyy, nznsczz)
    dqv = zeros(nxnscxx, nynscyy, nznsczz)

    return SgsTendencies(dn, dq, dqv)
end

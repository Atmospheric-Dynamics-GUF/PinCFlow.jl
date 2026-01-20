struct IceTendencies2{A <: AbstractArray{<:AbstractFloat, 3}}
    dn2::A
    dq2::A
    dqv2::A
end

function IceTendencies2(domain::Domain, ice_setup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    dn2 = zeros(nxx, nyy, nzz)
    dq2 = zeros(nxx, nyy, nzz)
    dqv2 = zeros(nxx, nyy, nzz)

    return IceTendencies2(dn, dq, dqv)
end

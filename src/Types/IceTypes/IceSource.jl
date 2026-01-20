struct IceSource{A <: AbstractArray{<:AbstractFloat, 3}}
    nsource :: A
    qsource :: A
    qvsource :: A
end

function IceSource(namelists::Namelists, domain::Domain)    
    (; ice_setup) = namelists.ice

    return IceSource(domain, ice_setup)
end

function IceSource(domain::Domain, ice_setup::NoIce)
    
    nsource = zeros(0, 0, 0)
    qsource = zeros(0, 0, 0)
    qvsource = zeros(0, 0, 0)

    return IceSource(nsource, qsource, qvsource)
end

function IceSource(domain::Domain, ice_setup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    nsource = zeros(nxx, nyy, nzz)
    qsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)

    return IceSource(nsource, qsource, qvsource)
end

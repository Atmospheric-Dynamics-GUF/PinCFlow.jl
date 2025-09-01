# module IceSource    

# using ..Types

# include("Ice.jl")
# include("IceAuxiliaries.jl")
# include("IcePredictands.jl")
# include("IceReconstructions.jl")
# include("IceTendencies.jl")
# include("IceFluxes.jl")

# export Ice, IceAuxiliaries, IcePredictands, IceReconstructions, IceTendencies, IceFluxes    
# end

struct IceSource{A <: AbstractArray{<:AbstractFloat, 3}}
    nsource :: A
    qsource :: A
    qvsource :: A
end

function IceSource(namelists::Namelists, domain::Domain)    
    (; icesetup) = namelists.ice

    return IceSource(domain, icesetup)
end

function IceSource(domain::Domain, icesetup::NoIce)
    
    nsource = zeros(0, 0, 0)
    qsource = zeros(0, 0, 0)
    qvsource = zeros(0, 0, 0)

    return IceSource(nsource, qsource, qvsource)
end

function IceSource(domain::Domain, icesetup::AbstractIce)
    (; nxx, nyy, nzz) = domain

    nsource = zeros(nxx, nyy, nzz)
    qsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)

    return IceSource(nsource, qsource, qvsource)
end

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

function IceSource(domain::Domain)
    (; nxx, nyy, nzz) = domain

    nsource = zeros(nxx, nyy, nzz)
    qsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)

    return IceSource(nsource, qsource, qvsource)
end

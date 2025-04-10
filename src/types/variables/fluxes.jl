struct Fluxes{A <: AbstractArray{<:AbstractFloat, 4}}
    phirho::A
    phirhop::A
    phiu::A
    phiv::A
    phiw::A
end

function Fluxes(domain::Domain)

    # Get parameters.
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw) = (zeros(nxx, nyy, nzz, 3) for i in 1:5)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw)
end

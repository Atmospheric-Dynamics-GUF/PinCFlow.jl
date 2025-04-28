struct Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
}
    phirho::A
    phirhop::A
    phiu::A
    phiv::A
    phiw::A
    phip::B
end

function Fluxes(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Fluxes(domain, model)
end

function Fluxes(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw) = (zeros(nxx, nyy, nzz, 3) for i in 1:5)
    phip = zeros(0, 0, 0, 0)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

function Fluxes(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw, phip) =
        (zeros(nxx, nyy, nzz, 3) for i in 1:6)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

"""
    Fluxes{A, B}

Storage for numerical fluxes in three spatial directions (x, y, z) on a 3D grid.

# Fields

  - `phirho::A`: Density flux arrays (nxx × nyy × nzz × 3)
  - `phirhop::A`: Density perturbation flux arrays (nxx × nyy × nzz × 3)
  - `phiu::A`: x-velocity flux arrays (nxx × nyy × nzz × 3)
  - `phiv::A`: y-velocity flux arrays (nxx × nyy × nzz × 3)
  - `phiw::A`: z-velocity flux arrays (nxx × nyy × nzz × 3)
  - `phip::B`: Pressure flux arrays (model-dependent dimensions)

The fourth dimension represents the three spatial directions for flux computation (x=1, y=2, z=3).
For incompressible models, `phip` is empty (0×0×0×0).
"""
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

"""
    Fluxes(namelists::Namelists, domain::Domain)

Construct `Fluxes` from configuration namelists and domain specification.
"""
function Fluxes(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Fluxes(domain, model)
end

"""
    Fluxes(domain::Domain, model::AbstractModel)

Construct `Fluxes` for incompressible models with empty pressure flux array.
"""
function Fluxes(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw) = (zeros(nxx, nyy, nzz, 3) for i in 1:5)
    phip = zeros(0, 0, 0, 0)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

"""
    Fluxes(domain::Domain, model::Compressible)

Construct `Fluxes` for compressible models including pressure flux arrays.
"""
function Fluxes(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw, phip) =
        (zeros(nxx, nyy, nzz, 3) for i in 1:6)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

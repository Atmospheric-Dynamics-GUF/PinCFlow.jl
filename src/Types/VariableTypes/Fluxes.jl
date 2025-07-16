"""
```julia
Fluxes{
    A <: AbstractArray{<:AbstractFloat, 4},
    B <: AbstractArray{<:AbstractFloat, 4},
}
```

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
```julia
Fluxes(namelists::Namelists, domain::Domain)
```

Construct `Fluxes` from configuration namelists and domain specification.

# Arguments

  - `namelists`: Configuration settings for the model
  - `domain`: Domain specification containing grid dimensions

# Returns

  - `::Fluxes`: `Fluxes` instance with zero-initialized arrays of the appropriate dimensions.
"""
function Fluxes(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Fluxes(domain, model)
end

"""
```julia
Fluxes(domain::Domain, model::AbstractModel)
```

Construct `Fluxes` for incompressible models with empty pressure flux array.

# Arguments

  - `domain`: Domain specification containing grid dimensions
  - `model`: Model type

# Returns

  - `::Fluxes`: `Fluxes` instance with zero-initialized arrays (`phip` has size `(0, 0, 0, 0)`).
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
```julia
Fluxes(domain::Domain, model::Compressible)
```

Construct `Fluxes` for compressible models including pressure flux arrays.

# Arguments

  - `domain`: Domain specification containing grid dimensions
  - `model`: Compressible model type

# Returns

  - `::Fluxes`: `Fluxes` instance with zero-initialized arrays.
"""
function Fluxes(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the fluxes.
    (phirho, phirhop, phiu, phiv, phiw, phip) =
        (zeros(nxx, nyy, nzz, 3) for i in 1:6)

    # Return a Fluxes instance.
    return Fluxes(phirho, phirhop, phiu, phiv, phiw, phip)
end

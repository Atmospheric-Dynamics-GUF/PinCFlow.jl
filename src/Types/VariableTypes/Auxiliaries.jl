"""
    Auxiliaries{A}

Container for auxiliary variables on a 3D grid.

# Fields

  - `phi::A`: Auxiliary scalar field (nxx × nyy × nzz)
"""
struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    phi::A
end

"""
    Auxiliaries(domain::Domain)

Construct `Auxiliaries` with zero-initialized auxiliary field based on domain dimensions.
"""
function Auxiliaries(domain::Domain)
    (; nxx, nyy, nzz) = domain

    # Initialize the auxiliaries.
    phi = zeros(nxx, nyy, nzz)

    # Return an Auxiliaries instance.
    return Auxiliaries(phi)
end

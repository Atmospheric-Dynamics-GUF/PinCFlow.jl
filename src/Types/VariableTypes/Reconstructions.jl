"""
```julia
Reconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
```

Arrays for the reconstructions of prognostic variables.

The first three dimensions represent physical space, the fourth dimension represents the direction in which the reconstruction was performed and the fifth dimension represents the two cell edges of the reconstruction.

# Fields

  - `rhotilde::A`: Reconstructed density.
  - `rhoptilde::A`: Reconstructed density fluctuations.
  - `utilde::A`: Reconstructed zonal momentum.
  - `vtilde::A`: Reconstructed meridional momentum.
  - `wtilde::A`: Reconstructed vertical momentum.
"""
struct Reconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    rhotilde::A
    rhoptilde::A
    utilde::A
    vtilde::A
    wtilde::A
end

"""
```julia
Reconstructions(domain::Domain)
```

Initialize arrays for the reconstruction of prognostic variables.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Returns

  - `::Reconstructions`: `Reconstructions` instance with zero-initialized arrays.
"""
function Reconstructions(domain::Domain)

    # Get parameters.
    (; nxx, nyy, nzz) = domain

    # Initialize the reconstructed variables.
    (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
        (zeros(nxx, nyy, nzz, 3, 2) for i in 1:5)

    # Return a Reconstructions instance.
    return Reconstructions(rhotilde, rhoptilde, utilde, vtilde, wtilde)
end

"""
    Reconstructions{A<:AbstractArray{<:AbstractFloat, 5}}

Storage for MUSCL-reconstructed variables at cell interfaces.

Holds left and right reconstructed states for all prognostic variables at cell faces
in each spatial direction.

# Type Parameters

  - `A<:AbstractArray{<:AbstractFloat, 5}`: 5D array type for reconstructed fields

# Fields

  - `rhotilde::A`: Reconstructed total density [ρ_L, ρ_R] at interfaces
  - `rhoptilde::A`: Reconstructed density perturbation [ρ'_L, ρ'_R] at interfaces
  - `utilde::A`: Reconstructed zonal velocity [u_L, u_R] at interfaces
  - `vtilde::A`: Reconstructed meridional velocity [v_L, v_R] at interfaces
  - `wtilde::A`: Reconstructed vertical velocity [w_L, w_R] at interfaces

# Array Structure

Each field has dimensions `(nxx, nyy, nzz, 3, 2)`:

  - Dimensions 1-3: Spatial grid points (i, j, k)
  - Dimension 4: Spatial direction (1=x, 2=y, 3=z)
  - Dimension 5: Interface side (1=left, 2=right)

# Constructor

    Reconstructions(domain::Domain)

Create reconstruction storage arrays sized for the computational domain.

# Usage Context

```julia
# Initialize reconstruction storage
recons = Reconstructions(domain)

# Apply MUSCL reconstruction
apply_3d_muscl!(state.variables.predictands, recons, state)

# Use reconstructed values in flux calculation
compute_fluxes!(state, recons)

# Access reconstructed values
u_left = recons.utilde[i, j, k, 1, 1]   # Left state, x-direction
u_right = recons.utilde[i, j, k, 1, 2]  # Right state, x-direction
```

# Physical Interpretation

  - Left/right states represent extrapolated values at cell interfaces
  - Enable higher-order spatial accuracy in flux calculations
  - Account for spatial gradients within computational cells
  - Support upwind-biased reconstruction for stability

See also: [`apply_3d_muscl!`](@ref), [`compute_fluxes!`](@ref), [`Variables`](@ref)
"""
struct Reconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
    rhotilde::A
    rhoptilde::A
    utilde::A
    vtilde::A
    wtilde::A
end

"""
    Reconstructions(domain::Domain)

Create MUSCL reconstruction storage for prognostic variables.

Allocates 5D arrays to store left and right reconstructed states at cell interfaces
for all prognostic fields in each spatial direction.

# Arguments

  - `domain::Domain`: Computational domain specification for array sizing

# Returns

  - `Reconstructions`: Initialized storage with zero-filled arrays

# Array Layout

  - Each field: `(nxx, nyy, nzz, 3, 2)` dimensions
  - Includes halo regions for boundary communications
  - Memory allocated for all spatial directions and interface sides
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

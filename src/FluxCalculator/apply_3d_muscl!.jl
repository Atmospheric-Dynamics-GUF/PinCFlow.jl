"""
```julia
apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limitertype::MCVariant,
)
```

Applies the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) reconstruction in three dimensions.

# Arguments

  - `phi::AbstractArray{<:AbstractFloat, 3}`: Input scalar values at cell centers.
  - `phitilde::AbstractArray{<:AbstractFloat, 5}`: Output reconstructed values at cell interfaces.
    The fourth dimension (1,2,3) represents directions (x,y,z), and the fifth dimension contains
    left/right interface values.
  - `nxx::Integer`: Size of the grid in x-direction, including ghost cells.
  - `nyy::Integer`: Size of the grid in y-direction, including ghost cells.
  - `nzz::Integer`: Size of the grid in z-direction, including ghost cells.
  - `limitertype::MCVariant`: Type of slope limiter to use.

# Details

This function applies one-dimensional MUSCL reconstruction sequentially in each direction
(x, y, and z). For each direction, it uses the `apply_1d_muscl!` function to perform
the reconstruction along the corresponding dimension.

# Notes

  - The function modifies `phitilde` in-place.
  - Reconstruction is applied only to interior cells (from 2 to n-1 in each dimension).
"""
function apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limitertype::MCVariant,
)

    # Reconstruct in x.
    for kz in 2:(nzz - 1), jy in 2:(nyy - 1)
        @views apply_1d_muscl!(phi[:, jy, kz], phitilde[:, jy, kz, 1, :], nxx)
    end

    # Reconstruct in y.
    for kz in 2:(nzz - 1), ix in 2:(nxx - 1)
        @views apply_1d_muscl!(phi[ix, :, kz], phitilde[ix, :, kz, 2, :], nyy)
    end

    # Reconstruct in z.
    for jy in 2:(nyy - 1), ix in 2:(nxx - 1)
        @views apply_1d_muscl!(phi[ix, jy, :], phitilde[ix, jy, :, 3, :], nzz)
    end

    # Return.
    return
end

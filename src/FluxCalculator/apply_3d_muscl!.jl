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

Applies the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) for reconstruction in three dimensions.

# Arguments

  - `phi`: Input array.
  - `phitilde`: Output array with reconstructed values. The fourth dimension represents the directions in which the input was reconstructed and the fifth dimension the reconstructions to the left and right.
  - `nxx`: Size of `phi` in ``\\widehat{x}``-direction.
  - `nyy`: Size of `phi` in ``\\widehat{y}``-direction.
  - `nzz`: Size of `phi` in ``\\widehat{z}``-direction.
  - `limitertype`: Type of flux limiter to use.

# See also

  - [`PinCFlow.FluxCalculator.apply_1d_muscl!`](@ref)
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

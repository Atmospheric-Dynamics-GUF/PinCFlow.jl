"""
```julia
apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limitertype::AbstractLimiter,
)
```

Apply the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) for reconstruction in three dimensions.

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
function apply_3d_muscl! end

function apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limitertype::AbstractLimiter,
)

    # Reconstruct in x.
    @ivy for k in 2:(nzz - 1), j in 2:(nyy - 1)
        apply_1d_muscl!(phi[:, j, k], phitilde[:, j, k, 1, :], nxx, limitertype)
    end

    # Reconstruct in y.
    @ivy for k in 2:(nzz - 1), i in 2:(nxx - 1)
        apply_1d_muscl!(phi[i, :, k], phitilde[i, :, k, 2, :], nyy, limitertype)
    end

    # Reconstruct in z.
    @ivy for j in 2:(nyy - 1), i in 2:(nxx - 1)
        apply_1d_muscl!(phi[i, j, :], phitilde[i, j, :, 3, :], nzz, limitertype)
    end

    return
end

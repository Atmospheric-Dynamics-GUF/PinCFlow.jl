"""
```julia
apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limiter_type::Val{:MCVariant},
)::Nothing
```

Apply the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) for reconstruction in three dimensions.

# Arguments

  - `phi`: Input array.

  - `phitilde`: Output array with reconstructed values. The fourth dimension represents the directions in which the input was reconstructed and the fifth dimension the reconstructions to the left and right.

  - `nxx`: Size of `phi` in ``\\hat{x}``-direction.

  - `nyy`: Size of `phi` in ``\\hat{y}``-direction.

  - `nzz`: Size of `phi` in ``\\hat{z}``-direction.

  - `limiter_type`: Type of flux limiter to use.

# See also

  - [`PinCFlow.FluxCalculator.apply_1d_muscl!`](@ref)
"""
function apply_3d_muscl! end

@ivy function apply_3d_muscl!(
    phi::AbstractArray{<:AbstractFloat, 3},
    phitilde::AbstractArray{<:AbstractFloat, 5},
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    limiter_type::Val{:MCVariant},
)::Nothing

    # Reconstruct in x.
    @fuse for k in 2:(nzz - 1), j in 2:(nyy - 1), i in 2:(nxx - 1)
        (left, right) = muscl(
            phi[i - 1, j, k],
            phi[i, j, k],
            phi[i + 1, j, k],
            limiter_type,
        )
        phitilde[i, j, k, 1, 1] = left
        phitilde[i, j, k, 1, 2] = right
    end

    # Reconstruct in y.
    @fuse for k in 2:(nzz - 1), j in 2:(nyy - 1), i in 2:(nxx - 1)
        (back, front) = muscl(
            phi[i, j - 1, k],
            phi[i, j, k],
            phi[i, j + 1, k],
            limiter_type,
        )
        phitilde[i, j, k, 2, 1] = back
        phitilde[i, j, k, 2, 2] = front
    end

    # Reconstruct in z.
    @fuse for k in 2:(nzz - 1), j in 2:(nyy - 1), i in 2:(nxx - 1)
        (down, up) = muscl(
            phi[i, j, k - 1],
            phi[i, j, k],
            phi[i, j, k + 1],
            limiter_type,
        )
        phitilde[i, j, k, 3, 1] = down
        phitilde[i, j, k, 3, 2] = up
    end

    nothing
end

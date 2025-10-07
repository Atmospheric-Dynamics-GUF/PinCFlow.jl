"""
```julia
apply_1d_muscl!(
    phi::AbstractVector{<:AbstractFloat},
    phitilde::AbstractMatrix{<:AbstractFloat},
    phisize::Integer,
    limiter_type::MCVariant,
)
```

Apply the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) for reconstruction in one dimension.

The reconstruction to the left is given by

```math
{\\widetilde{\\phi}}^\\mathrm{L} = \\begin{cases}
    \\phi & \\mathrm{if} \\quad \\phi = \\phi_{i - 1} \\quad \\mathrm{or} \\quad \\phi = \\phi_{i + 1},\\\\
    \\phi - \\frac{1}{2} \\eta \\left(\\frac{\\phi_{i + 1} - \\phi}{\\phi - \\phi_{i - 1}}\\right) \\left(\\phi - \\phi_{i - 1}\\right) & \\mathrm{else}
\\end{cases}
```

and that to the right is given by

```math
{\\widetilde{\\phi}}^\\mathrm{R} = \\begin{cases}
    \\phi & \\mathrm{if} \\quad \\phi = \\phi_{i - 1} \\quad \\mathrm{or} \\quad \\phi = \\phi_{i + 1},\\\\
    \\phi + \\frac{1}{2} \\eta \\left(\\frac{\\phi - \\phi_{i - 1}}{\\phi_{i + 1} - \\phi}\\right) \\left(\\phi_{i + 1} - \\phi\\right) & \\mathrm{else},
\\end{cases}
```

where

```math
\\eta \\left(\\xi\\right) = \\max \\left[0, \\min \\left(2 \\xi, \\frac{2 + \\xi}{3}, 2\\right)\\right]
```

is the monotonized-centered variant limiter.

# Arguments

  - `phi`: Input vector.

  - `phitilde`: Output matrix with reconstructed values. The two columns of `phitilde` contain the reconstructions to the left and right. No reconstruction is computed for the first and last row of `phitilde`.

  - `phisize`: Length of the input vector `phi`.

  - `limiter_type`: Type of flux limiter to use.
"""
function apply_1d_muscl! end

function apply_1d_muscl!(
    phi::AbstractVector{<:AbstractFloat},
    phitilde::AbstractMatrix{<:AbstractFloat},
    phisize::Integer,
    limiter_type::MCVariant,
)

    # Initialize phitilde.
    phitilde .= 1000.0

    # Reconstruct.
    @ivy for i in 2:(phisize - 1)
        deltal = phi[i] - phi[i - 1]
        deltar = phi[i + 1] - phi[i]

        if deltar == 0.0
            phitilde[i, 2] = phi[i]
            phitilde[i, 1] = phi[i]
        else
            if deltal == 0.0
                theta = deltal / deltar
                s = (2.0 + theta) / 3.0
                sigmal = max(0.0, min(2 * theta, s, 2.0))

                phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
                phitilde[i, 1] = phi[i]
            else
                theta = deltal / deltar

                s = (2.0 + theta) / 3.0
                sigmal = max(0.0, min(2 * theta, s, 2.0))

                s = (2.0 + 1.0 / theta) / 3.0
                sigmar = max(0.0, min(2 / theta, s, 2.0))

                phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
                phitilde[i, 1] = phi[i] - 0.5 * sigmar * deltal
            end
        end
    end

    return
end

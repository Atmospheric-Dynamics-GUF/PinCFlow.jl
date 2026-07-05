"""
```julia
muscl(
    phil::AbstractFloat,
    phic::AbstractFloat,
    phir::AbstractFloat,
    limiter_type::Val{:MCVariant},
)::NTuple{2, <:AbstractFloat}
```

Apply the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) for reconstruction in one dimension.

The reconstruction to the left is given by

```math
{\\tilde{\\phi}}^\\mathrm{L} = \\begin{cases}
    \\phi & \\mathrm{if} \\quad \\phi = \\phi_{i - 1} \\quad \\mathrm{or} \\quad \\phi = \\phi_{i + 1},\\\\
    \\phi - \\frac{1}{2} \\eta \\left(\\frac{\\phi_{i + 1} - \\phi}{\\phi - \\phi_{i - 1}}\\right) \\left(\\phi - \\phi_{i - 1}\\right) & \\mathrm{else}
\\end{cases}
```

and that to the right is given by

```math
{\\tilde{\\phi}}^\\mathrm{R} = \\begin{cases}
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

  - `phil`: Value to the left.

  - `phic`: Value in the center.

  - `phir`: Value to the right.

  - `limiter_type`: Type of flux limiter to use.
"""
function muscl end

function muscl(
    phil::AbstractFloat,
    phic::AbstractFloat,
    phir::AbstractFloat,
    limiter_type::Val{:MCVariant},
)::NTuple{2, <:AbstractFloat}
    deltal = phic - phil
    deltar = phir - phic

    if deltar == 0.0
        return (phic, phic)
    else
        if deltal == 0.0
            theta = deltal / deltar
            s = (2.0 + theta) / 3.0
            sigmal = max(0.0, min(2 * theta, s, 2.0))

            return (phic, phic + 0.5 * sigmal * deltar)
        else
            theta = deltal / deltar

            s = (2.0 + theta) / 3.0
            sigmal = max(0.0, min(2 * theta, s, 2.0))

            s = (2.0 + 1.0 / theta) / 3.0
            sigmar = max(0.0, min(2 / theta, s, 2.0))

            return (phic - 0.5 * sigmar * deltal, phic + 0.5 * sigmal * deltar)
        end
    end
end

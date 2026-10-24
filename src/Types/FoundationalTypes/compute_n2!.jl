"""
```julia
compute_n2!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    thetabar::AbstractArray{<:AbstractFloat, 3},
    n2::AbstractArray{<:AbstractFloat, 3},
)::Nothing
```

Compute the buoyancy frequency ``N^2 \\left(z\\right)`` from the potential temperature.

The squared buoyancy frequency is given by

```math
\\begin{align*}
N^2 & = \\frac{g}{\\bar{\\theta}} \\frac{\\bar{\\theta}_{k + 1} - \\bar{\\theta}_{k - 1}}{2 J \\Delta \\hat{z}} \\;.
\\end{align*}
```

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `thetabar`: Potential-temperature background.

  - `n2`: Squared buoyancy frequency.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_boundaries_of_field!`](@ref)
"""
function compute_n2! end

@ivy function compute_n2!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    thetabar::AbstractArray{<:AbstractFloat, 3},
    n2::AbstractArray{<:AbstractFloat, 3},
)::Nothing
    (; z_size, nbz) = namelists.domain
    (; g_ndim) = constants
    (; nz, ko, k0, k1, nxx, nyy) = domain
    (; jac, dz) = grid

    # Compute the squared buoyancy frequency.
    @fuse n2 = 0.0
    @fuse for k in k0:k1, j in 1:nyy, i in 1:nxx
        n2[i, j, k] =
            g_ndim / thetabar[i, j, k] / jac[i, j, k] *
            0.5 *
            (thetabar[i, j, k + 1] - thetabar[i, j, k - 1]) / dz
    end

    # Compute the squared buoyancy frequency at the boundaries.
    set_vertical_boundaries_of_field!(n2, namelists, domain, +)
    if ko == 0
        @fuse for k in 1:nbz, j in 1:nyy, i in 1:nxx
            n2[i, j, k] =
                g_ndim / thetabar[i, j, k0 - 1] / jac[i, j, k0 - 1] *
                (thetabar[i, j, k0] - thetabar[i, j, k0 - 1]) / dz
        end
    end
    if ko + nz == z_size
        @fuse for k in 1:nbz, j in 1:nyy, i in 1:nxx
            n2[i, j, k1 + k] =
                g_ndim / thetabar[i, j, k1 + 1] / jac[i, j, k1 + 1] *
                (thetabar[i, j, k1 + 1] - thetabar[i, j, k1]) / dz
        end
    end

    nothing
end

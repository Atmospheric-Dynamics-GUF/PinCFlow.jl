"""
```julia 
compute_n2!(
)
```

Compute the buoyancy frequency ``N^2 \\left(z\\right)`` from the potential temperature.

The squared buoyancy frequency is given by

```math 
\\begin{align*}
N^2 & = \\frac{g}{\\overline{\\theta}} \\frac{\\overline{\\theta}_{k + 1} - \\overline{\\theta}_{k - 1}}{2 J \\Delta \\widehat{z}} \\;.
\\end{align*}
```

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.
  
# See also

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_boundaries_of_field!`](@ref)
"""
function compute_n2! end

function compute_n2!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    thetabar::AbstractArray{<:AbstractFloat, 3},
    n2::AbstractArray{<:AbstractFloat, 3},
)
    (; nbz) = namelists.domain
    (; g_ndim) = constants
    (; zz_size, nzz, ko, k0, k1) = domain
    (; jac, dz) = grid

    # Compute the squared buoyancy frequency.
    n2 .= 0.0
    @ivy for k in k0:k1
        n2[:, :, k] .=
            g_ndim ./ thetabar[:, :, k] ./ jac[:, :, k] .* 0.5 .*
            (thetabar[:, :, k + 1] .- thetabar[:, :, k - 1]) ./ dz
    end

    # Compute the squared buoyancy frequency at the boundaries.
    set_vertical_boundaries_of_field!(n2, namelists, domain, +)
    @ivy if ko == 0
        for k in 1:nbz
            n2[:, :, k] .=
                g_ndim ./ thetabar[:, :, k0 - 1] ./ jac[:, :, k0 - 1] .*
                (thetabar[:, :, k0] .- thetabar[:, :, k0 - 1]) ./ dz
        end
    end
    @ivy if ko + nzz == zz_size
        for k in 1:nbz
            n2[:, :, k1 + k] .=
                g_ndim ./ thetabar[:, :, k1 + 1] ./ jac[:, :, k1 + 1] .*
                (thetabar[:, :, k1 + 1] .- thetabar[:, :, k1]) ./ dz
        end
    end

    return
end
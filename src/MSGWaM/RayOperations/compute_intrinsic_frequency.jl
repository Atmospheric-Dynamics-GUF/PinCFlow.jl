"""
```julia
compute_intrinsic_frequency(state::State, indices::NTuple{4, <:Integer})
```

Compute the intrinsic frequency of the ray volume specified by `indices`.

Calculates the intrinsic frequency from the dispersion relation

```math
\\widehat{\\omega}_\\alpha = \\sigma \\sqrt{\\frac{N^2 \\left(k_\\alpha^2 + l_\\alpha^2\\right) + f^2 m_\\alpha^2}{\\left|\\boldsymbol{k}_\\alpha\\right|^2}},
```

where ``\\boldsymbol{k}_\\alpha = \\left(k_\\alpha, l_\\alpha, m_\\alpha\\right)^\\mathrm{T}`` is the ray volume's wavevector, ``N^2`` is the squared buoyancy frequency interpolated to the ray volume's vertical position, ``f`` is the Coriolis parameter and ``\\sigma`` is the frequency branch.

# Arguments

  - `state`: Model state

  - `indices`: Indices of the ray volume of interest.

# Returns

  - `::AbstractFloat`: Intrinsic frequency.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)
"""
function compute_intrinsic_frequency end

function compute_intrinsic_frequency(
    state::State,
    indices::NTuple{4, <:Integer},
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; rays) = state.wkb

    zr = rays.z[indices...]
    kr = rays.k[indices...]
    lr = rays.l[indices...]
    mr = rays.m[indices...]
    khr = sqrt(kr^2 + lr^2)

    n2r = interpolate_stratification(zr, state, N2())
    fc = coriolis_frequency * tref

    return branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)
end

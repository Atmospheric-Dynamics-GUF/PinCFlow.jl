"""
```julia
compute_intrinsic_frequency(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
```

Return the intrinsic frequency of the ray volume specified by ``\\left(r, i, j, k\\right)``.

The intrinsic frequency is calculated from the dispersion relation

```math
\\widehat{\\omega}_r = \\sigma \\sqrt{\\frac{N_r^2 \\left(k_r^2 + l_r^2\\right) + f^2 m_r^2}{\\left|\\boldsymbol{k}_r\\right|^2}},
```

where ``N_r^2`` is the squared buoyancy frequency interpolated to the ray volume's vertical position and ``\\sigma`` is the frequency branch.

# Arguments

  - `state`: Model state

  - `r`: Ray-volume index.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)
"""
function compute_intrinsic_frequency end

function compute_intrinsic_frequency(
    state::State,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; tref) = state.constants
    (; rays) = state.wkb

    @ivy zr = rays.z[r, i, j, k]
    @ivy kr = rays.k[r, i, j, k]
    @ivy lr = rays.l[r, i, j, k]
    @ivy mr = rays.m[r, i, j, k]
    khr = sqrt(kr^2 + lr^2)

    n2r = interpolate_stratification(zr, state, N2())
    fc = coriolis_frequency * tref

    return branch * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)
end

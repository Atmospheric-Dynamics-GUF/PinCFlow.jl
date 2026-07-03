"""
```julia
compute_blocked_layer!(
    state::State,
    deltah::AbstractFloat,
    n2h::AbstractFloat,
    uh::AbstractFloat,
    vh::AbstractFloat,
    i::Integer,
    j::Integer,
)::AbstractFloat
```

Update the blocked-layer depth and return the ratio between the effective and total mountain-wave amplitude.

The blocked-layer depth is given by

```math
\\Delta z_\\mathrm{B} = 2 \\Delta h \\max \\left(0, \\frac{\\mathrm{Lo} - C}{\\mathrm{Lo}}\\right),
```

where ``\\mathrm{Lo} = N_h \\Delta h \\left|\\boldsymbol{k}_h\\right| / \\left|\\boldsymbol{k}_h \\cdot \\boldsymbol{u}_h\\right|`` represents the Long number (with ``\\boldsymbol{k}_h`` computed by `compute_slope`) and ``C`` the Long-number threshold (`state.namelists.wkb.long_threshold`). The corresponding ratio between effective and total mountain-wave amplitude is given by

```math
r \\left(\\Delta z_\\mathrm{B}\\right) = 1 - B \\frac{\\Delta z_\\mathrm{B}}{2 \\Delta h}
```

with ``B`` regulating how much of the blocked layer contributes to the amplitude reduction (`state.namelists.wkb.reduction_coefficient`).

# Arguments

  - `state`: Model state.

  - `deltah`: Elevation difference between the local background orography and the summits of the true local orography.

  - `n2h`: Squared buoyancy frequency averaged between ``h_{\\mathrm{b}}`` and ``h_{\\mathrm{b}} + \\Delta h``.

  - `uh`: Zonal wind averaged between ``h_{\\mathrm{b}}`` and ``h_{\\mathrm{b}} + \\Delta h``.

  - `vh`: Meridional wind averaged between ``h_{\\mathrm{b}}`` and ``h_{\\mathrm{b}} + \\Delta h``.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.

# See also

  - [`PinCFlow.MSGWaM.BlockedLayer.compute_slope`](@ref)

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
function compute_blocked_layer! end

@ivy function compute_blocked_layer!(
    state::State,
    deltah::AbstractFloat,
    n2h::AbstractFloat,
    uh::AbstractFloat,
    vh::AbstractFloat,
    i::Integer,
    j::Integer,
)::AbstractFloat
    (; blocking, long_threshold, reduction_coefficient) = state.namelists.wkb
    (; k0) = state.domain
    (; zctilde) = state.grid
    (; deltazb) = state.wkb

    if blocking && deltah > 0
        kh = compute_slope(state, deltah, i, j)
        deltazb[i, j] =
            2 *
            deltah *
            max(
                0,
                1 -
                long_threshold / sqrt(n2h) / deltah / sqrt(kh[1]^2 + kh[2]^2) *
                abs(kh[1] * uh + kh[2] * vh),
            )
        ratio = 1 - reduction_coefficient * deltazb[i, j] / 2 / deltah
    elseif blocking
        ratio = 1
        deltazb[i, j] = 0
    else
        ratio = 1
    end

    return ratio
end

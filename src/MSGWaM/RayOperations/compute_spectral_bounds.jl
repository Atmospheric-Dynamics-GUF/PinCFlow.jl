"""
```julia
compute_spectral_bounds(
    wavenumbers::AbstractVector{<:AbstractFloat},
)::NTuple{4, <:AbstractFloat}
```

Compute the minima and maxima of positive and negative entries in `wavenumbers` and return them (in the order positive minimum, positive maximum, negative minimum, negative maximum).

This method is used by [`PinCFlow.MSGWaM.RayUpdate.merge_rays!`](@ref) to create spectral bins.

# Arguments

  - `wavenumbers`: Vector of wavenumbers in the considered spectral dimension.
"""
function compute_spectral_bounds end

function compute_spectral_bounds(
    wavenumbers::AbstractVector{<:AbstractFloat},
)::NTuple{4, <:AbstractFloat}

    # Compute minima and maxima.
    (krminp, krmaxp) =
        extrema(kr for kr in wavenumbers if kr > 0; init = (0.0, 0.0))
    (krminn, krmaxn) =
        extrema(kr for kr in wavenumbers if kr < 0; init = (0.0, 0.0))

    # Adjust bounds that are zero.
    if krminp == krmaxp == krminn == krmaxn == 0
        krminp = 1.0
        krmaxp = 2.0
        krminn = 1.0
        krmaxn = 2.0
    elseif krminp == krmaxp == 0
        krminp = krminn
        krmaxp = krmaxn
    elseif krminn == krmaxn == 0
        krminn = krminp
        krmaxn = krmaxp
    end

    # Prevent zero-width intervals.
    if krminn == krmaxn
        krminn /= 2
        krmaxn *= 2
    end
    if krminp == krmaxp
        krminp /= 2
        krmaxp *= 2
    end

    (krminp, krmaxp, krminn, krmaxn)
end

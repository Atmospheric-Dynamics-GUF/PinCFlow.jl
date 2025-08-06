"""
```julia
compute_merge_index(
    wnr::AbstractFloat,
    wnr_min_p::AbstractFloat,
    wnr_max_p::AbstractFloat,
    wnr_min_n::AbstractFloat,
    wnr_max_n::AbstractFloat,
    dwnr_mrg_p::AbstractFloat,
    dwnr_mrg_n::AbstractFloat,
    nray::Integer,
)
```

Computes the index of the wavenumber `wnr` on a 1D spectral grid specified by the other arguments.

This method is used by [`PinCFlow.MSGWaM.RayUpdate.merge_rays!`](@ref) to sort ray volumes into spectral bins.

# Arguments

  - `wnr`: Wavenumber value.
  - `wnr_min_p`: Minimum positive wavenumber.
  - `wnr_max_p`: Maximum positive wavenumber.
  - `wnr_min_n`: Minimum negative wavenumber.
  - `wnr_max_n`: Maximum negative wavenumber.
  - `dwnr_mrg_p`: Logarithmic spacing of discrete positive wavenumbers.
  - `dwnr_mrg_n`: Logarithmic spacing of discrete negative wavenumbers.
  - `nray`: Number of spectral grid points.

# Returns

  - `::Integer`: Position on the 1D spectral grid.
"""
function compute_merge_index end

function compute_merge_index(
    wnr::AbstractFloat,
    wnr_min_p::AbstractFloat,
    wnr_max_p::AbstractFloat,
    wnr_min_n::AbstractFloat,
    wnr_max_n::AbstractFloat,
    dwnr_mrg_p::AbstractFloat,
    dwnr_mrg_n::AbstractFloat,
    nray::Integer,
)
    if wnr < 0
        if abs(log(-wnr / wnr_max_n) / dwnr_mrg_n) < 1
            iray = div(nray, 2) - 1
        else
            iray = round(Int, log(-wnr / wnr_min_n) / dwnr_mrg_n) + 1
        end
    elseif wnr == 0
        iray = div(nray, 2)
    else
        if abs(log(wnr / wnr_max_p) / dwnr_mrg_p) < 1
            iray = nray - 1
        else
            iray =
                round(Int, log(wnr / wnr_min_p) / dwnr_mrg_p) + div(nray, 2) + 1
        end
    end

    return iray
end

"""
    compute_merge_index(wnr::AbstractFloat, wnr_min_p::AbstractFloat, wnr_max_p::AbstractFloat, wnr_min_n::AbstractFloat, wnr_max_n::AbstractFloat, dwnr_mrg_p::AbstractFloat, dwnr_mrg_n::AbstractFloat, nray::Integer)

# Arguments

  - `wnr::AbstractFloat`: Wavenumber value
  - `wnr_min_p::AbstractFloat`: Minimum positive wavenumber
  - `wnr_max_p::AbstractFloat`: Maximum positive wavenumber
  - `wnr_min_n::AbstractFloat`: Minimum negative wavenumber
  - `wnr_max_n::AbstractFloat`: Maximum negative wavenumber
  - `dwnr_mrg_p::AbstractFloat`: Positive wavenumber merge spacing
  - `dwnr_mrg_n::AbstractFloat`: Negative wavenumber merge spacing
  - `nray::Integer`: Number of ray volumes
"""
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

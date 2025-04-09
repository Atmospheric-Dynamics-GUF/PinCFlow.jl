function compute_merge_index(
  wnr::AbstractFloat,
  wnr_min_p::AbstractFloat,
  wnr_max_p::AbstractFloat,
  wnr_min_n::AbstracFloat,
  wnr_max_n::AbstractFloat,
  dwnr_mrg_p::AbstractFloat,
  dwnr_mrg_n::AbstractFloat,
  nray::AbstractFloat,
)
  if wnr < 0
    if abs(log(-wnr / wnr_max_n) / dwnr_mrg_n) < 1
      iray = nray / 2 - 1
    else
      iray = round(Int, log(-wnr / wnr_min_n) / dwnr_mrg_n) + 1
    end
  elseif wnr == 0
    iray = nray / 2
  else
    if abs(log(wnr / wnr_max_p) / dwnr_mrg_p) < 1
      iray = nray - 1
    else
      iray = round(Int, log[wnr / wnr_min_p] / dwnr_mrg_p) + nray / 2 + 1
    end
  end

  return iray
end

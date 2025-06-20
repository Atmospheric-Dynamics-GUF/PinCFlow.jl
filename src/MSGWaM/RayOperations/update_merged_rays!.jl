"""
    update_merged_rays!(merge_mode::AbstractMergeMode, merged_rays::AbstractVector{<:MergedRays}, jray::Integer, xr::AbstractFloat, dxr::AbstractFloat, yr::AbstractFloat, dyr::AbstractFloat, zr::AbstractFloat, dzr::AbstractFloat, kr::AbstractFloat, dkr::AbstractFloat, lr::AbstractFloat, dlr::AbstractFloat, mr::AbstractFloat, dmr::AbstractFloat, fxk::AbstractFloat, fyl::AbstractFloat, fzm::AbstractFloat, nr::AbstractFloat, omegar::AbstractFloat)

Update merged ray volume data structure during ray merging process.

Accumulates ray volume information into spectral bins, expanding spatial and
spectral extents as needed and integrating wave action densities.

# Arguments

  - `merge_mode::AbstractMergeMode`: Merging strategy (affects wave action integration)
  - `merged_rays::AbstractVector{<:MergedRays}`: Array of merged ray bins to update
  - `jray::Integer`: Bin index for this ray volume
  - `xr, yr, zr::AbstractFloat`: Ray physical position coordinates
  - `dxr, dyr, dzr::AbstractFloat`: Ray physical extent in each direction
  - `kr, lr, mr::AbstractFloat`: Ray spectral position coordinates
  - `dkr, dlr, dmr::AbstractFloat`: Ray spectral extent in each direction
  - `fxk, fyl, fzm::AbstractFloat`: Phase space volume factors
  - `nr::AbstractFloat`: Wave action density
  - `omegar::AbstractFloat`: Intrinsic frequency

# Merging Process

 1. **First Ray in Bin**: Initialize extent boundaries from ray edges
 2. **Subsequent Rays**: Expand boundaries to encompass all rays in bin
 3. **Wave Action Integration**: Accumulate wave action using merge mode-specific formula

# Spatial Extent Expansion

For each direction (x, y, z):

  - `min_boundary = min(current_min, ray_center - ray_extent/2)`
  - `max_boundary = max(current_max, ray_center + ray_extent/2)`

# Spectral Extent Expansion

For each wavenumber (k, l, m):

  - `min_boundary = min(current_min, k_center - k_extent/2)`
  - `max_boundary = max(current_max, k_center + k_extent/2)`

# Wave Action Accumulation

Uses `compute_wave_action_integral` with merge mode to properly weight
the contribution based on:

  - Raw wave action density `nr`
  - Frequency factor `omegar`
  - Phase space volume factors `fxk`, `fyl`, `fzm`

# Applications

Called during the ray merging process for each ray being combined into
a spectral bin, building up the merged ray characteristics that will
be used to create the final merged ray volumes.
"""
function update_merged_rays!(
    merge_mode::AbstractMergeMode,
    merged_rays::AbstractVector{<:MergedRays},
    jray::Integer,
    xr::AbstractFloat,
    dxr::AbstractFloat,
    yr::AbstractFloat,
    dyr::AbstractFloat,
    zr::AbstractFloat,
    dzr::AbstractFloat,
    kr::AbstractFloat,
    dkr::AbstractFloat,
    lr::AbstractFloat,
    dlr::AbstractFloat,
    mr::AbstractFloat,
    dmr::AbstractFloat,
    fxk::AbstractFloat,
    fyl::AbstractFloat,
    fzm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    if merged_rays[jray].nr[] == 0
        merged_rays[jray].xr .= [xr - dxr / 2, xr + dxr / 2]
        merged_rays[jray].yr .= [yr - dyr / 2, yr + dyr / 2]
        merged_rays[jray].zr .= [zr - dzr / 2, zr + dzr / 2]
        merged_rays[jray].kr .= [kr - dkr / 2, kr + dkr / 2]
        merged_rays[jray].lr .= [lr - dlr / 2, lr + dlr / 2]
        merged_rays[jray].mr .= [mr - dmr / 2, mr + dmr / 2]
    else
        merged_rays[jray].xr .= [
            min(merged_rays[jray].xr[1], xr - dxr / 2),
            max(merged_rays[jray].xr[2], xr + dxr / 2),
        ]
        merged_rays[jray].yr .= [
            min(merged_rays[jray].yr[1], yr - dyr / 2),
            max(merged_rays[jray].yr[2], yr + dyr / 2),
        ]
        merged_rays[jray].zr .= [
            min(merged_rays[jray].zr[1], zr - dzr / 2),
            max(merged_rays[jray].zr[2], zr + dzr / 2),
        ]
        merged_rays[jray].kr .= [
            min(merged_rays[jray].kr[1], kr - dkr / 2),
            max(merged_rays[jray].kr[2], kr + dkr / 2),
        ]
        merged_rays[jray].lr .= [
            min(merged_rays[jray].lr[1], lr - dlr / 2),
            max(merged_rays[jray].lr[2], lr + dlr / 2),
        ]
        merged_rays[jray].mr .= [
            min(merged_rays[jray].mr[1], mr - dmr / 2),
            max(merged_rays[jray].mr[2], mr + dmr / 2),
        ]
    end

    merged_rays[jray].nr[] +=
        compute_wave_action_integral(merge_mode, nr, omegar, fxk, fyl, fzm)

    return
end

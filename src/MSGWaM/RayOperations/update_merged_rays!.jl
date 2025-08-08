"""
```julia
update_merged_rays!(
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
```

Update the fields of `merged_rays[jray]` such that they contain the outermost bounds and total wave action/energy of all contributing ray volumes.

This method is used to compute the properties of merged ray volumes. It is called for every old ray volume that contributes to the new, merged volume and updates the outermost bounds in physical and spectral space, as well as the total wave action/energy, accordingly.

# Arguments

- `merge_mode`: Merging strategy.
- `merged_rays`: Array of merged ray volumes.
- `jray`: Index of the merged ray volume to update.
- `xr`: Position of the old ray volume in ``x``.
- `dxr`: Extent of the old ray volume in ``x``.
- `yr`: Position of the old ray volume in ``y``.
- `dyr`: Extent of the old ray volume in ``y``.
- `zr`: Position of the old ray volume in ``z``.
- `dzr`: Extent of the old ray volume in ``z``.
- `kr`: Position of the old ray volume in ``k``.
- `dkr`: Extent of the old ray volume in ``k``.
- `lr`: Position of the old ray volume in ``l``.
- `dlr`: Extent of the old ray volume in ``k``.
- `mr`: Position of the old ray volume in ``m``.
- `dmr`: Extent of the old ray volume in ``m``.
- `fxk`: Phase-space factor of the old ray volume in ``x``-``k`` subspace.
- `fyl`: Phase-space factor of the old ray volume in ``y``-``l`` subspace.
- `fzm`: Phase-space factor of the old ray volume in ``z``-``m`` subspace.
- `nr`: Phase-space wave-action density of the old ray volume.
- `omegar`: Intrinsic frequency of the old ray volume.

# See also

- [`PinCFlow.MSGWaM.RayOperations.compute_wave_action_integral`](@ref)
"""
function update_merged_rays! end

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

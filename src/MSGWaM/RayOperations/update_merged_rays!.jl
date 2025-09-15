"""
```julia
update_merged_rays!(
    merge_mode::AbstractMergeMode,
    merged_rays::MergedRays,
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

Update the fields of `merged_rays` at `jray` such that they contain the outermost bounds and total wave action/energy of all contributing ray volumes.

This method is used to compute the properties of merged ray volumes. It is called for every old ray volume that contributes to the new, merged volume and updates the outermost bounds in physical and spectral space, as well as the total wave action/energy, accordingly.

# Arguments

  - `merge_mode`: Merging strategy.

  - `merged_rays`: Properties of merged ray volumes.

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
    merged_rays::MergedRays,
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
    @ivy if merged_rays.nr[jray] == 0
        for (i, o) in ((1, -), (2, +))
            merged_rays.xr[i, jray] = o(xr, dxr / 2)
            merged_rays.yr[i, jray] = o(yr, dyr / 2)
            merged_rays.zr[i, jray] = o(zr, dzr / 2)
            merged_rays.kr[i, jray] = o(kr, dkr / 2)
            merged_rays.lr[i, jray] = o(lr, dlr / 2)
            merged_rays.mr[i, jray] = o(mr, dmr / 2)
        end
    else
        for (i, e, o) in ((1, min, -), (2, max, +))
            merged_rays.xr[i, jray] = e(merged_rays.xr[i, jray], o(xr, dxr / 2))
            merged_rays.yr[i, jray] = e(merged_rays.yr[i, jray], o(yr, dyr / 2))
            merged_rays.zr[i, jray] = e(merged_rays.zr[i, jray], o(zr, dzr / 2))
            merged_rays.kr[i, jray] = e(merged_rays.kr[i, jray], o(kr, dkr / 2))
            merged_rays.lr[i, jray] = e(merged_rays.lr[i, jray], o(lr, dlr / 2))
            merged_rays.mr[i, jray] = e(merged_rays.mr[i, jray], o(mr, dmr / 2))
        end
    end

    merged_rays.nr[jray] +=
        compute_wave_action_integral(merge_mode, nr, omegar, fxk, fyl, fzm)

    return
end

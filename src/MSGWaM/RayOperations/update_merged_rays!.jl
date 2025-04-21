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

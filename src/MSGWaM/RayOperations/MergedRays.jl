struct MergedRays{A <: SVector{2, <:AbstractFloat}, B <: AbstractFloat}
    xr::A
    yr::A
    zr::A
    kr::A
    lr::A
    mr::A
    nr::B
end

function MergedRays()
    return MergedRays([SVector{2}(0.0, 0.0) for i in 1:7]...)
end

function MergedRays(
    merge_mode::AbstractMergeMode,
    self::MergedRays,
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
    axk::AbstractFloat,
    ayl::AbstractFloat,
    azm::AbstractFloat,
    nr::AbstractFloat,
    omegar::AbstractFloat,
)
    if self.nr == 0
        return MergedRays(
            SVector{2}(xr - dxr / 2, xr + dxr / 2),
            SVector{2}(yr - dyr / 2, yr + dyr / 2),
            SVector{2}(zr - dzr / 2, zr + dzr / 2),
            SVector{2}(kr - dkr / 2, kr + dkr / 2),
            SVector{2}(lr - dlr / 2, lr + dlr / 2),
            SVector{2}(mr - dmr / 2, mr + dmr / 2),
            merge_wave_action(merge_mode, axk, ayl, azm, nr, omegar),
        )
    else
        return MergedRays(
            SVector{2}(
                min(self.xr[1], xr - dxr / 2),
                max(self.xr[2], xr + dxr / 2),
            ),
            SVector{2}(
                min(self.yr[1], yr - dyr / 2),
                max(self.yr[2], yr + dyr / 2),
            ),
            SVector{2}(
                min(self.zr[1], zr - dzr / 2),
                max(self.zr[2], zr + dzr / 2),
            ),
            SVector{2}(
                min(self.kr[1], kr - dkr / 2),
                max(self.kr[2], kr + dkr / 2),
            ),
            SVector{2}(
                min(self.lr[1], lr - dlr / 2),
                max(self.lr[2], lr + dlr / 2),
            ),
            SVector{2}(
                min(self.mr[1], mr - dmr / 2),
                max(self.mr[2], mr + dmr / 2),
            ),
            merge_wave_action(merge_mode, self, axk, ayl, azm, nr, omegar),
        )
    end
end

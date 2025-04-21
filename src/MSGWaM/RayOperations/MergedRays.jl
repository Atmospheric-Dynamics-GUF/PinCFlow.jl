struct MergedRays{
    A <: AbstractVector{<:AbstractFloat},
    B <: Ref{<:AbstractFloat},
}
    xr::A
    yr::A
    zr::A
    kr::A
    lr::A
    mr::A
    nr::B
end

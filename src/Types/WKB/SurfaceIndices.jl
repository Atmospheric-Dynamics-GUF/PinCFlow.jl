struct SurfaceIndices{
    A <: AbstractArray{<:Integer, 3},
    B <: AbstractVector{<:Integer},
}
    ir_sfc::A
    ix2_sfc::B
    jy2_sfc::B
    kz2_sfc::B
    ik_sfc::B
    jl_sfc::B
    km_sfc::B
    iwm_sfc::B
end

function SurfaceIndices(n_sfc::Integer, nxx::Integer, nyy::Integer)
    return SurfaceIndices(
        zeros(Int, n_sfc, nxx, nyy),
        [zeros(Int, n_sfc) for i in 1:7]...,
    )
end

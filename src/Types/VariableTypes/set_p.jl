function set_p(
    model::AbstractModel,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)
    return zeros(0, 0, 0)
end

function set_p(
    model::Compressible,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
    pstrattfc::AbstractArray{<:AbstractFloat, 3},
)
    return copy(pstrattfc)
end
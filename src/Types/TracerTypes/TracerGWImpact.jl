struct TracerGWImpact{A <: AbstractArray{<:AbstractFloat, 3}}
    uchi::A
    vchi::A
    wchi::A
    dchidt::A
end

function TracerGWImpact(nxi::Integer, nyi::Integer, nzi::Integer)
    return TracerGWImpact([zeros(nxi, nyi, nzi) for i in 1:4]...)
end
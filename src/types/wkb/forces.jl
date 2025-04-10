struct Forces{A <: AbstractArray{<:AbstractFloat, 3}}
    u::A
    v::A
    w::A
end

function Forces(nxx::Integer, nyy::Integer, nzz::Integer)
    return Forces([zeros(nxx, nyy, nzz) for i in 1:3]...)
end
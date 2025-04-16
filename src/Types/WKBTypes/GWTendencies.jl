struct GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dudt::A
    dvdt::A
    dthetadt::A
end

function GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end

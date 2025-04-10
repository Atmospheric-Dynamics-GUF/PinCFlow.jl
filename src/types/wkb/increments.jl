struct Increments{A <: AbstractArray{<:AbstractFloat, 4}}
    dxray::A
    dyray::A
    dzray::A
    dkray::A
    dlray::A
    dmray::A
    ddxray::A
    ddyray::A
    ddzray::A
end

function Increments(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
    return Increments([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:9]...)
end

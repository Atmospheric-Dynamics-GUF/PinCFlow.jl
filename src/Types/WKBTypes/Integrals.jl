struct Integrals{A <: AbstractArray{<:AbstractFloat, 3}}
    uu::A
    uv::A
    uw::A
    vv::A
    vw::A
    etx::A
    ety::A
    utheta::A
    vtheta::A
    e::A
    dudt::A
    dvdt::A
    dthetadt::A
end

function Integrals(nxx::Integer, nyy::Integer, nzz::Integer)
    return Integrals([zeros(nxx, nyy, nzz) for i in 1:13]...)
end

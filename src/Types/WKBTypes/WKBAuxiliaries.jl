struct WKBAuxiliaries{
    A <: AbstractArray{<:ComplexF64, 3},
    B <: AbstractArray{ComplexF64, 2},
}
    uhatold::A
    vhatold::A
    whatold::A
    bhatold::A
    pihatold::A
    mat::B
end

function WKBAuxiliaries(
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::WKBAuxiliaries
    return WKBAuxiliaries(
        [zeros(ComplexF64, nxx, nyy, nzz) for i in 1:5]...,
        zeros(ComplexF64, 5, 5),
    )
end

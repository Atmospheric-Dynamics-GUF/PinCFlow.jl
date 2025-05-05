struct BicGStab{
    A <: AbstractMatrix{<:AbstractFloat},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    r_vm::A
    p::B
    r0::B
    rold::B
    r::B
    s::B
    t::B
    v::B
    matvec::B
    v_pc::B
end

function BicGStab(namelists::Namelists, domain::Domain)
    (; sizex, sizey) = namelists.domain
    (; nx, ny, nz) = domain
    return BicGStab(zeros(sizex, sizey), [zeros(nx, ny, nz) for i in 1:9]...)
end

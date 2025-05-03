struct Preconditioner{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractMatrix{<:AbstractFloat},
}
    s_pc::A
    q_pc::A
    p_pc::B
    s_pc_bc::B
    q_pc_bc::B
end

function Preconditioner(domain::Domain)

    # Get all necessary fields.
    (; nx, ny, nz) = domain

    # Return a Preconditioner instance.
    return Preconditioner(
        [zeros(nx, ny, nz) for i in 1:2]...,
        [zeros(nx, ny) for i in 1:3]...,
    )
end

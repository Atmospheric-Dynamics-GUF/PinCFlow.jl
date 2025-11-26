function thomas_algorithm! end

function thomas_algorithm!(
    state::State,
    ath::AbstractArray{<:AbstractFloat, 3},
    bth::AbstractArray{<:AbstractFloat, 3},
    cth::AbstractArray{<:AbstractFloat, 3},
    fth::AbstractArray{<:AbstractFloat, 3},
    qth::AbstractArray{<:AbstractFloat, 3},
    pth::AbstractMatrix{<:AbstractFloat},
    fth_bc::AbstractMatrix{<:AbstractFloat},
    qth_bc::AbstractMatrix{<:AbstractFloat},
)
    (; comm, nz, nz, ko, up, down) = state.domain
    (; z_size) = state.namelists.domain

    if ko == 0
        qth[:, :, 1] .= .-cth[:, :, 1] ./ bth[:, :, 1]
        fth[:, :, 1] .= fth[:, :, 1] ./ bth[:, :, 1]
    else
        MPI.Recv!(qth_bc, comm; source = down, tag = 1)
        MPI.Recv!(fth_bc, comm; source = down, tag = 2)

        pth .= 1.0 ./ (bth[:, :, 1] .+ ath[:, :, 1] .* qth_bc)
        qth[:, :, 1] .= .-cth[:, :, 1] .* pth
        fth[:, :, 1] .= (fth[:, :, 1] .- ath[:, :, 1] .* fth_bc) .* pth
    end

    for k in 2:nz
        pth .= 1.0 ./ (bth[:, :, k] .+ ath[:, :, k] .* qth[:, :, k - 1])
        qth[:, :, k] .= .-cth[:, :, k] .* pth
        fth[:, :, k] .=
            (fth[:, :, k] .- ath[:, :, k] .* fth[:, :, k - 1]) .* pth
    end

    if ko + nz != z_size
        qth_bc .= qth[:, :, nz]
        fth_bc .= fth[:, :, nz]

        MPI.Send(qth_bc, comm; dest = up, tag = 1)
        MPI.Send(fth_bc, comm; dest = up, tag = 2)

        MPI.Recv!(fth_bc, comm; source = up)

        fth[:, :, nz] .+= qth[:, :, nz] .* fth_bc
    end

    for k in (nz - 1):-1:1
        fth[:, :, k] .+= qth[:, :, k] .* fth[:, :, k + 1]
    end

    if ko != 0
        fth_bc .= fth[:, :, 1]

        MPI.Send(fth_bc, comm; dest = down)
    end

    return
end

"""
```julia 
thomas_algorithm!(state::State)
```

Solves a tridiagonal system in ``\\hat{z}``-direction using the Thomas tridiagonal matrix algorithm (see [Durran, 2010](https://doi.org/10.1007/978-1-4419-6412-0)) . Since the Thomas algorithm consists of an upward elimination sweep and a downward pass, this method performs sequential one-way MPI communication if the domain is parallelized in the vertical.

The system is defined as:

```math 
a_k \\phi_{k-1} + b_k\\phi_k + c_k\\phi_{k+1} = f_k\\;.
```

The result is stored in `state.variables.auxiliaries.fth`.

# Arguments 

  - `state`: Model state.
"""
function thomas_algorithm! end

function thomas_algorithm!(state::State)
    (; comm, nz, ko, up, down) = state.domain
    (; z_size) = state.namelists.domain
    (; ath, bth, cth, fth, qth, pth, fth_bc, qth_bc) =
        state.variables.auxiliaries

    @ivy if ko == 0
        qth[:, :, 1] .= .-cth[:, :, 1] ./ bth[:, :, 1]
        fth[:, :, 1] .= fth[:, :, 1] ./ bth[:, :, 1]
    else
        MPI.Recv!(qth_bc, comm; source = down, tag = 1)
        MPI.Recv!(fth_bc, comm; source = down, tag = 2)

        pth .= 1.0 ./ (bth[:, :, 1] .+ ath[:, :, 1] .* qth_bc)
        qth[:, :, 1] .= .-cth[:, :, 1] .* pth
        fth[:, :, 1] .= (fth[:, :, 1] .- ath[:, :, 1] .* fth_bc) .* pth
    end

    @ivy for k in 2:nz
        pth .= 1.0 ./ (bth[:, :, k] .+ ath[:, :, k] .* qth[:, :, k - 1])
        qth[:, :, k] .= .-cth[:, :, k] .* pth
        fth[:, :, k] .=
            (fth[:, :, k] .- ath[:, :, k] .* fth[:, :, k - 1]) .* pth
    end

    @ivy if ko + nz != z_size
        qth_bc .= qth[:, :, nz]
        fth_bc .= fth[:, :, nz]

        MPI.Send(qth_bc, comm; dest = up, tag = 1)
        MPI.Send(fth_bc, comm; dest = up, tag = 2)

        MPI.Recv!(fth_bc, comm; source = up)

        fth[:, :, nz] .+= qth[:, :, nz] .* fth_bc
    end

    @ivy for k in (nz - 1):-1:1
        fth[:, :, k] .+= qth[:, :, k] .* fth[:, :, k + 1]
    end

    @ivy if ko != 0
        fth_bc .= fth[:, :, 1]

        MPI.Send(fth_bc, comm; dest = down)
    end

    return
end

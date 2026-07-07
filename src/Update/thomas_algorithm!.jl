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

@ivy function thomas_algorithm!(state::State)
    (; comm, nx, ny, nz, ko, up, down) = state.domain
    (; z_size) = state.namelists.domain
    (; ath, bth, cth, fth, qth, fth_bc, qth_bc) = state.variables.auxiliaries

    if ko == 0
        @share for j in 1:ny, i in 1:nx
            qth[i, j, 1] = -cth[i, j, 1] / bth[i, j, 1]
            fth[i, j, 1] = fth[i, j, 1] / bth[i, j, 1]
        end
    else
        MPI.Recv!(qth_bc, comm; source = down, tag = 1)
        MPI.Recv!(fth_bc, comm; source = down, tag = 2)

        @share for j in 1:ny, i in 1:nx
            pth = 1.0 / (bth[i, j, 1] + ath[i, j, 1] * qth_bc[i, j])
            qth[i, j, 1] = -cth[i, j, 1] * pth
            fth[i, j, 1] = (fth[i, j, 1] - ath[i, j, 1] * fth_bc[i, j]) * pth
        end
    end

    for k in 2:nz
        @share for j in 1:ny, i in 1:nx
            pth = 1.0 / (bth[i, j, k] + ath[i, j, k] * qth[i, j, k - 1])
            qth[i, j, k] = -cth[i, j, k] * pth
            fth[i, j, k] =
                (fth[i, j, k] - ath[i, j, k] * fth[i, j, k - 1]) * pth
        end
    end

    if ko + nz != z_size
        @share for j in 1:ny, i in 1:nx
            qth_bc[i, j] = qth[i, j, nz]
            fth_bc[i, j] = fth[i, j, nz]
        end

        MPI.Send(qth_bc, comm; dest = up, tag = 1)
        MPI.Send(fth_bc, comm; dest = up, tag = 2)

        MPI.Recv!(fth_bc, comm; source = up)

        @share fth[:, :, nz] += qth[:, :, nz] * fth_bc
    end

    for k in (nz - 1):-1:1
        @share fth[:, :, k] += qth[:, :, k] * fth[:, :, k + 1]
    end

    if ko != 0
        @share fth_bc = fth[:, :, 1]

        MPI.Send(fth_bc, comm; dest = down)
    end

    return
end

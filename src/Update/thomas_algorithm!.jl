"""
```julia 
thomas_algorithm!(state::State)
```

Solves a tridiagonal system in ``\\hat{z}``-direction by dispatching to the appropriate vertical boundary condition method.

```julia
thomas_algorithm!(
    state::State,
    vertical_boundary_condition::Val{:Periodic},
)
```

Solves a tridiagonal system in ``\\hat{z}``-direction for periodic vertical boundary conditions using the Thomas tridiagonal matrix algorithm (see [Durran, 2010](https://doi.org/10.1007/978-1-4419-6412-0)) . Since the Thomas algorithm consists of an upward elimination sweep and a downward pass, this method performs sequential one-way MPI communication if the domain is parallelized in the vertical.

The system is defined as:

```math
a_k \\phi_{k-1} + b_k\\phi_k + c_k\\phi_{k+1} = f_k\\;.
```

The result is stored in `state.variables.auxiliaries.fth`.

```julia
thomas_algorithm!(
    state::State,
    vertical_boundary_condition::Val{:SolidWall},
)
```

Solves a tridiagonal system in ``\\hat{z}``-direction for solid wall vertical boundary conditions using the Thomas tridiagonal matrix algorithm (see [Durran, 2010](https://doi.org/10.1007/978-1-4419-6412-0)) . Since the Thomas algorithm consists of an upward elimination sweep and a downward pass, this method performs sequential one-way MPI communication if the domain is parallelized in the vertical.

The system is defined as for the method above with the result stored in `state.variables.auxiliaries.fth`.

# Arguments

  - `state`: Model state.
"""
function thomas_algorithm! end

function thomas_algorithm!(state::State)
    (; vertical_boundary_condition) = state.namelists.domain

    @dispatch_vertical_boundary_condition thomas_algorithm!(
        state,
        Val(vertical_boundary_condition),
    )

    return
end

@ivy function thomas_algorithm!(
    state::State,
    vertical_boundary_condition::Val{:Periodic},
)
    (; comm, nz, ko, up, down) = state.domain
    (; z_size) = state.namelists.domain
    (; ath, bth, cth, fth, qth, sth, pth, fth_bc, qth_bc, sth_bc, fthnz) =
        state.variables.auxiliaries

    fthnz .= fth[:, :, nz]

    if ko == 0
        qth[:, :, 1] .= .-cth[:, :, 1] ./ bth[:, :, 1]
        fth[:, :, 1] .= fth[:, :, 1] ./ bth[:, :, 1]
        sth[:, :, 1] .= .-ath[:, :, 1] ./ bth[:, :, 1]
    else
        MPI.Recv!(qth_bc, comm; source = down, tag = 1)
        MPI.Recv!(fth_bc, comm; source = down, tag = 2)
        MPI.Recv!(sth_bc, comm; source = down, tag = 3)

        pth .= 1.0 ./ (bth[:, :, 1] .+ ath[:, :, 1] .* qth_bc)
        qth[:, :, 1] .= .-cth[:, :, 1] .* pth
        fth[:, :, 1] .= (fth[:, :, 1] .- ath[:, :, 1] .* fth_bc) .* pth
        sth[:, :, 1] .= .-ath[:, :, 1] .* sth_bc .* pth
    end

    for k in 2:nz
        pth .= 1.0 ./ (bth[:, :, k] .+ ath[:, :, k] .* qth[:, :, k - 1])
        qth[:, :, k] .= .-cth[:, :, k] .* pth
        fth[:, :, k] .=
            (fth[:, :, k] .- ath[:, :, k] .* fth[:, :, k - 1]) .* pth
        sth[:, :, k] .= .-ath[:, :, k] .* sth[:, :, k - 1] .* pth
    end

    if ko + nz != z_size
        qth_bc .= qth[:, :, nz]
        fth_bc .= fth[:, :, nz]
        sth_bc .= sth[:, :, nz]

        MPI.Send(qth_bc, comm; dest = up, tag = 1)
        MPI.Send(fth_bc, comm; dest = up, tag = 2)
        MPI.Send(sth_bc, comm; dest = up, tag = 3)
    end

    if ko + nz == z_size
        sth[:, :, nz] .= 1.0
        qth[:, :, nz] .= 0.0
    else
        MPI.Recv!(sth_bc, comm; source = up, tag = 1)
        MPI.Recv!(qth_bc, comm; source = up, tag = 2)

        sth[:, :, nz] .+= qth[:, :, nz] .* sth_bc
        qth[:, :, nz] .= fth[:, :, nz] .+ qth[:, :, nz] .* qth_bc
    end

    for k in (nz - 1):-1:1
        sth[:, :, k] .+= qth[:, :, k] .* sth[:, :, k + 1]
        qth[:, :, k] .= fth[:, :, k] .+ qth[:, :, k] .* qth[:, :, k + 1]
    end

    if ko != 0
        sth_bc .= sth[:, :, 1]
        qth_bc .= qth[:, :, 1]

        MPI.Send(sth_bc, comm; dest = down, tag = 1)
        MPI.Send(qth_bc, comm; dest = down, tag = 2)
    end

    if ko + nz == z_size
        fth[:, :, nz] .=
            (
                fthnz .- cth[:, :, nz] .* qth[:, :, 1] .-
                ath[:, :, nz] .* qth[:, :, nz - 1]
            ) ./ (
                cth[:, :, nz] .* sth[:, :, 1] .+
                ath[:, :, nz] .* sth[:, :, nz - 1] .+ bth[:, :, nz]
            )
        
        fth_bc .= fth[:, :, nz]

        MPI.Bcast!(fth_bc, comm; root = MPI.Comm_rank(comm))
    end

    for k in 1:(nz - 1)
        fth[:, :, k] .= fth_bc .* sth[:, :, k] .+ qth[:, :, k]
    end

    return
end

@ivy function thomas_algorithm!(
    state::State,
    vertical_boundary_condition::Val{:SolidWall},
)
    (; comm, nz, ko, up, down) = state.domain
    (; z_size) = state.namelists.domain
    (; ath, bth, cth, fth, qth, pth, fth_bc, qth_bc) =
        state.variables.auxiliaries

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

"""
```julia
apply_preconditioner!(
    sin::AbstractArray{<:AbstractFloat, 3},
    sout::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    poisson::Poisson,
)
```

Apply line relaxation preconditioner for the Poisson equation.

This preconditioner uses alternating direction implicit (ADI) method with vertical
line relaxation to accelerate convergence of the BiCGStab solver. It treats the
vertical direction implicitly while horizontal coupling is handled explicitly.

# Arguments

  - `sin`: Input residual field
  - `sout`: Preconditioned output field
  - `namelists`: Contains preconditioner parameters (dtau, maxiteradi)
  - `domain`: MPI domain decomposition info for vertical communication
  - `grid`: Grid spacing for pseudo-time step calculation
  - `poisson`: Operator coefficients and preconditioner workspace

# Algorithm

 1. Set pseudo-time step based on horizontal grid spacing
 2. Iterate ADI relaxation sweeps
 3. Perform tridiagonal solves in vertical direction
 4. Handle MPI communication for domain boundaries
 5. Apply upward and downward elimination sweeps

# See also

  - [`PinCFlow.PoissonSolver.apply_operator!`](@ref)
"""
function apply_preconditioner!(
    sin::AbstractArray{<:AbstractFloat, 3},
    sout::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    grid::Grid,
    poisson::Poisson,
)
    (; dtau, maxiteradi) = namelists.poisson
    (; comm, sizezz, nz, nzz, ko, down, up) = domain
    (; dx, dy) = grid
    (; au_b, ac_b, ad_b) = poisson.tensor
    (; s_pc, q_pc, p_pc, s_pc_bc, q_pc_bc) = poisson.preconditioner

    # Initialize auxiliary fields.
    s_pc .= 0.0
    q_pc .= 0.0
    p_pc .= 0.0

    # Set pseudo-time step.
    deta = dtau / (2 * (1 / dx^2 + 1 / dy^2))

    # Iterate.
    for niter in 1:maxiteradi
        apply_operator!(s_pc, q_pc, Horizontal(), namelists, domain, poisson)
        s_pc .+= deta .* (q_pc .- sin)

        # Set the lower boundary.
        if ko == 0
            @views q_pc[:, :, 1] .=
                deta .* au_b[:, :, 1] ./ (1 .- deta .* ac_b[:, :, 1])
            @views s_pc[:, :, 1] ./= 1 .- deta .* ac_b[:, :, 1]
        else
            MPI.Recv!(q_pc_bc, comm; source = down, tag = 1)
            MPI.Recv!(s_pc_bc, comm; source = down, tag = 2)

            @views p_pc .=
                1 ./
                (1 .- deta .* ac_b[:, :, 1] .- deta .* ad_b[:, :, 1] .* q_pc_bc)
            @views q_pc[:, :, 1] .= deta .* au_b[:, :, 1] .* p_pc
            @views s_pc[:, :, 1] .=
                (s_pc[:, :, 1] .+ deta .* ad_b[:, :, 1] .* s_pc_bc) .* p_pc
        end

        # Perform upward sweep.
        for k in 2:nz
            @views p_pc .=
                1 ./ (
                    1 .- deta .* ac_b[:, :, k] .-
                    deta .* ad_b[:, :, k] .* q_pc[:, :, k - 1]
                )
            @views q_pc[:, :, k] .= deta .* au_b[:, :, k] .* p_pc
            @views s_pc[:, :, k] .=
                (s_pc[:, :, k] .+ deta .* ad_b[:, :, k] .* s_pc[:, :, k - 1]) .*
                p_pc
        end

        # Communicate the upper boundary and set it for the downward sweep.
        if ko + nzz != sizezz
            @views q_pc_bc .= q_pc[:, :, nz]
            @views s_pc_bc .= s_pc[:, :, nz]

            MPI.Send(q_pc_bc, comm; dest = up, tag = 1)
            MPI.Send(s_pc_bc, comm; dest = up, tag = 2)

            MPI.Recv!(s_pc_bc, comm; source = up)

            @views s_pc[:, :, nz] .+= q_pc[:, :, nz] .* s_pc_bc
        end

        # Perform downward sweep.
        for k in (nz - 1):-1:1
            @views s_pc[:, :, k] .+= q_pc[:, :, k] .* s_pc[:, :, k + 1]
        end

        # Communicate the lower boundary.
        if ko != 0
            @views s_pc_bc .= s_pc[:, :, 1]

            MPI.Send(s_pc_bc, comm; dest = down)
        end
    end

    sout .= s_pc

    return
end

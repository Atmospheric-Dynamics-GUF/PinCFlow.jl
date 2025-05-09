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
        @views if ko == 0
            q_pc[:, :, 1] .=
                deta .* au_b[:, :, 1] ./ (1 .- deta .* ac_b[:, :, 1])
            s_pc[:, :, 1] ./= 1 .- deta .* ac_b[:, :, 1]
        else
            MPI.Recv!(q_pc_bc, comm; source = down, tag = 1)
            MPI.Recv!(s_pc_bc, comm; source = down, tag = 2)

            p_pc .=
                1 ./
                (1 .- deta .* ac_b[:, :, 1] .- deta .* ad_b[:, :, 1] .* q_pc_bc)
            q_pc[:, :, 1] .= deta .* au_b[:, :, 1] .* p_pc
            s_pc[:, :, 1] .=
                (s_pc[:, :, 1] .+ deta .* ad_b[:, :, 1] .* s_pc_bc) .* p_pc
        end

        # Perform upward sweep.
        @views for k in 2:nz
            p_pc .=
                1 ./ (
                    1 .- deta .* ac_b[:, :, k] .-
                    deta .* ad_b[:, :, k] .* q_pc[:, :, k - 1]
                )
            q_pc[:, :, k] .= deta .* au_b[:, :, k] .* p_pc
            s_pc[:, :, k] .=
                (s_pc[:, :, k] .+ deta .* ad_b[:, :, k] .* s_pc[:, :, k - 1]) .*
                p_pc
        end

        # Communicate the upper boundary and set it for the downward sweep.
        @views if ko + nzz != sizezz
            q_pc_bc .= q_pc[:, :, nz]
            s_pc_bc .= s_pc[:, :, nz]

            MPI.Send(q_pc_bc, comm; dest = up, tag = 1)
            MPI.Send(s_pc_bc, comm; dest = up, tag = 2)

            MPI.Recv!(s_pc_bc, comm; source = up)

            s_pc[:, :, nz] .+= q_pc[:, :, nz] .* s_pc_bc
        end

        # Perform downward sweep.
        @views for k in (nz - 1):-1:1
            s_pc[:, :, k] .+= q_pc[:, :, k] .* s_pc[:, :, k + 1]
        end

        # Communicate the lower boundary.
        @views if ko != 0
            s_pc_bc .= s_pc[:, :, 1]

            MPI.Send(s_pc_bc, comm; dest = down)
        end
    end

    # Set final result.
    sout .= s_pc

    # Return.
    return
end

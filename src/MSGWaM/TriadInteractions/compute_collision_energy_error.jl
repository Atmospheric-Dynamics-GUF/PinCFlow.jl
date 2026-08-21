function compute_collision_energy_error end

function compute_collision_energy_error(
    state::State,
    ::Triad2D;
    print_diagnostic::Bool = true,
)
    (; master, comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; x_size) = state.namelists.domain
    (; dx, dy, dz, jac) = state.grid
    (; n2) = state.atmosphere
    (; spec_tend) = state
    (; col_int) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid

    local_positive_transfer = 0.0
    local_negative_transfer = 0.0
    local_absolute_transfer = 0.0

    @ivy for kk in k0:k1, jj in j0:j1, ii in i0:i1

        # Physical-cell volume.
        physical_cell_volume = dx * dy * jac[ii, jj, kk] * dz

        # Local buoyancy frequency.
        nn = sqrt(n2[ii, jj, kk])

        for mi in eachindex(m), kpi in eachindex(kp)

            st = col_int[ii, jj, kk, kpi, mi]

            if iszero(st)
                continue
            end

            # Use the same normalized intrinsic frequency as the
            # kinetic-equation implementation.
            omega = nn * compute_omega_hat(kp[kpi], m[mi])

            # Horizontal spectral measure:
            #
            # x_size == 1 : discrete k sum, measure = 1
            # x_size > 1  : continuous k integral, measure = Δk
            spectral_measure = if x_size == 1
                delm[mi]
            else
                delkp[kpi] * delm[mi]
            end

            energy_transfer =
                omega * st *
                spectral_measure *
                physical_cell_volume

            if !isfinite(energy_transfer)
                error(
                    "Nonfinite collision-energy tendency detected at ",
                    "(ii,jj,kk,kpi,mi) = ",
                    (ii, jj, kk, kpi, mi),
                    ": ",
                    energy_transfer,
                )
            end

            if energy_transfer > 0.0
                local_positive_transfer += energy_transfer
            elseif energy_transfer < 0.0
                local_negative_transfer += energy_transfer
            end

            local_absolute_transfer += abs(energy_transfer)
        end
    end

    local_diagnostic = [
        local_positive_transfer,
        local_negative_transfer,
        local_absolute_transfer,
    ]

    global_diagnostic = MPI.Allreduce(local_diagnostic, +, comm)

    positive_transfer = global_diagnostic[1]
    negative_transfer = global_diagnostic[2]
    absolute_transfer = global_diagnostic[3]

    net_transfer = positive_transfer + negative_transfer

    energy_error =
        absolute_transfer > 0.0 ?
        net_transfer / absolute_transfer :
        0.0

    if master && print_diagnostic
        println("")
        println("Collision-integral energy-conservation diagnostic:")
        println("  Positive energy transfer       = ", positive_transfer)
        println("  Negative energy transfer       = ", negative_transfer)
        println("  Net energy tendency            = ", net_transfer)
        println("  Absolute energy transfer       = ", absolute_transfer)
        println("  Relative energy error          = ", energy_error)
        println("  Relative energy error (%)      = ", 100.0 * energy_error)
        println("")
    end

    return energy_error
end
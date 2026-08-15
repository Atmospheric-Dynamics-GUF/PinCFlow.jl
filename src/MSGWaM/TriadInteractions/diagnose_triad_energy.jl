function diagnose_triad_energy end

function diagnose_triad_energy(
    state::State,
    ::Triad2D,
)
    (; master, comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; x_size, y_size) = state.namelists.domain

    (; spec_tend) = state.wkb
    (; wavespectrum, was_pred) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid

    (; dx, dy, dz, jac) = state.grid
    (; n2) = state.atmosphere

    local_energy_before = 0.0
    local_energy_after = 0.0

    #----------------------------------------------------------
    # Eulerian wave energy before and after the triad update:
    #
    # E = ∫ ω̂(k,m) N(k,m) dk dm dV
    #
    # with hydrostatic intrinsic frequency
    #
    # ω̂ = N k / |m|.
    #----------------------------------------------------------

    for kk in k0:k1, jj in j0:j1, ii in i0:i1

        n_local = sqrt(n2[ii, jj, kk])

        dx_cell = x_size > 1 ? dx : 1.0
        dy_cell = y_size > 1 ? dy : 1.0

        physical_cell_volume =
            dx_cell * dy_cell * jac[ii, jj, kk] * dz

        for mi in eachindex(m), kpi in eachindex(kp)

            omega_hat =
                n_local * kp[kpi] / abs(m[mi])

            spectral_cell_width =
                abs(delkp[kpi] * delm[mi])

            energy_weight =
                omega_hat *
                spectral_cell_width *
                physical_cell_volume

            local_energy_before +=
                was_pred[ii, jj, kk, kpi, mi] *
                energy_weight

            local_energy_after +=
                wavespectrum[ii, jj, kk, kpi, mi] *
                energy_weight
        end
    end

    #----------------------------------------------------------
    # Global energies over all MPI subdomains.
    #----------------------------------------------------------

    energy_before =
        MPI.Allreduce(local_energy_before, +, comm)

    energy_after =
        MPI.Allreduce(local_energy_after, +, comm)

    energy_change =
        energy_after - energy_before

    relative_energy_change =
        energy_before > 0.0 ?
        energy_change / energy_before :
        0.0

    if master
        println("")
        println("Eulerian triad-energy diagnostic:")
        println("  Energy before interaction = ", energy_before)
        println("  Energy after interaction  = ", energy_after)
        println("  Energy change             = ", energy_change)
        println("  Relative energy change    = ", relative_energy_change)
        println("")
    end

    return (
        energy_before = energy_before,
        energy_after = energy_after,
        energy_change = energy_change,
        relative_energy_change = relative_energy_change,
    )
end
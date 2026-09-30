function compute_action_ref! end

function compute_action_ref!(
    state::State;
    support_tol::AbstractFloat = 0.0,
    diagonal_connectivity::Bool = false,
    verify::Bool = true,
)
    (; master, comm, i0, i1, j0, j1, k0, k1, ko) = state.domain
    (; x_size) = state.namelists.domain
    (; source_mode) = state.namelists.wkb

    (; spec_tend) = state
    (; wavespectrum, action_ref) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; lref) = state.constants

    nkp = length(kp)
    nm = length(m)

    local_peak_action = zeros(Float64, nkp, nm)

    @ivy for mi in eachindex(m), kpi in eachindex(kp)

        spectral_cell_measure =
            x_size == 1 ?
            abs(delm[mi]) :
            abs(delkp[kpi] * delm[mi])

        peak_action = 0.0

        if source_mode isa ContinuousSpectralSource

            # Only the bottom MPI rank owns the prescribed
            # continuous-source reservoir at k0 - 1.
            if ko == 0
                kk = k0 - 1

                for jj in j0:j1, ii in i0:i1
                    was = wavespectrum[ii, jj, kk, kpi, mi]

                    if isfinite(was) && was > 0.0
                        spectral_cell_action =
                            was * spectral_cell_measure

                        if spectral_cell_action > peak_action
                            peak_action = spectral_cell_action
                        end
                    end
                end
            end

        else

            # Ordinary initialization: identify modes from
            # the physical domain only.
            for kk in k0:k1, jj in j0:j1, ii in i0:i1
                was = wavespectrum[ii, jj, kk, kpi, mi]

                if isfinite(was) && was > 0.0
                    spectral_cell_action =
                        was * spectral_cell_measure

                    if spectral_cell_action > peak_action
                        peak_action = spectral_cell_action
                    end
                end
            end
        end

        local_peak_action[kpi, mi] = peak_action
    end

    global_peak_action =
        MPI.Allreduce(local_peak_action, max, comm)

    mode_information =
        get_initial_mode_information(
            global_peak_action;
            support_tol = support_tol,
            diagonal_connectivity = diagonal_connectivity,
        )

    if isempty(mode_information)
        error(
            "Error in compute_action_ref!: no initialized ",
            "spectral mode was identified. ",
            "support_tol = ",
            support_tol,
        )
    end

    if :chi_parent in state.namelists.output.output_variables
        initialize_chi_parent_tracker!(state, mode_information)
    end

    action_ref_value = Inf

    for mode in mode_information
        action_ref_value =
            min(action_ref_value, mode.peak_action)
    end

    action_ref[] = action_ref_value

    if verify && master
        print_initial_mode_information(
            mode_information,
            kp,
            m,
            lref,
            action_ref[],
            support_tol,
            diagonal_connectivity,
        )
    end

    return nothing
end
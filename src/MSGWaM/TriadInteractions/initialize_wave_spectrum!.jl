function initialize_wave_spectrum! end

function initialize_wave_spectrum!(state::State)
    (; wkb_mode) = state.namelists.wkb
    (; triad_mode) = state.namelists.triad
    initialize_wave_spectrum!(state, wkb_mode, triad_mode)
    return
end


function initialize_wave_spectrum!(state::State, 
    wkb_mode::Union{NoWKB, MultiColumn, SingleColumn, SteadyState}, 
    triad_mode::NoTriad)
    return
end

function initialize_wave_spectrum!(
    state::State,
    wkb_mode::Union{MultiColumn, SingleColumn, SteadyState},
    triad_mode::Union{Triad2D, Triad3DIso},
)
    (; master) = state.domain
    (; i0, i1, j0, j1, k0, k1) = state.domain

    (; nthreads_triad, compute_dephasing_time, action_rel_tol) =
        state.namelists.triad

    (; spec_tend) = state
    (; wavespectrum, nl_time_scale, dephasing_time, prev_dt, action_ref) =
        spec_tend

    (; kp, m, delkp, delm) = spec_tend.spec_grid

    if master
        println(repeat("-", 80))
        println("")
        println("Triad interaction module activated with the model: ", triad_mode)
        println("Number of threads per rank: (nthreads_triad) = ", nthreads_triad)
        println("")
        println(repeat("-", 80))
    end

    # Project the initialized rays onto the wave-spectrum grid.
    get_wave_spectrum!(state)
    
    # Compute the reference action once from the initialized
    # spectral modes.
    compute_action_ref!(state; support_tol = 0.0, diagonal_connectivity = false, verify = true)

    if !isfinite(action_ref[]) || action_ref[] <= 0.0
        error("Invalid action_ref = ", action_ref[])
    end

    # Minimum significant spectral-cell action.
    action_floor = action_rel_tol * action_ref[]

    # Reset the timescale arrays before calculating the initial
    # nonlinear and dephasing times.
    nl_time_scale .= Inf
    dephasing_time .= Inf
    
    if compute_dephasing_time && triad_mode isa Triad3DIso
        error("Dephasing-time calculation is currently implemented only for Triad2D.")
    end

    @ivy for kk in k0:k1, jj in j0:j1, ii in i0:i1

        #------------------------------------------------------
        # Check whether this physical cell contains at least
        # one materially active spectral cell.
        #------------------------------------------------------

        max_cell_action = 0.0

        for mi in eachindex(m), kpi in eachindex(kp)
            spectral_cell_width = delkp[kpi] * delm[mi]
            cell_action = wavespectrum[ii, jj, kk, kpi, mi] * spectral_cell_width

            if cell_action > max_cell_action
                max_cell_action = cell_action
            end
        end

        if max_cell_action <= action_floor
            continue
        end

        #------------------------------------------------------
        # Calculate the initial collision integral and nonlinear
        # timescale for this physical cell.
        #------------------------------------------------------

        compute_scattering_integral!(state, ii, jj, kk, triad_mode)

        tau_nl = get_nl_time_scale(spec_tend, ii, jj, kk, action_rel_tol)
        nl_time_scale[ii, jj, kk] = tau_nl

        #------------------------------------------------------
        # Calculate the initial dephasing timescale using the
        # same initial wave spectrum and collision integral.
        #------------------------------------------------------

        if compute_dephasing_time && triad_mode isa Triad2D
            tau_pl = get_dephasing_time(state, ii, jj, kk, tau_nl, triad_mode)
            dephasing_time[ii, jj, kk] = tau_pl
        end
    end
    
    # Disable timestep-growth restriction for the first step.
    prev_dt[] = Inf

    return nothing
end
function update_wave_spectrum! end

function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat,
    triad_mode::Union{Triad2D}
    )
    (; time_scheme) = state.namelists.triad
    update_wave_spectrum!(state, ii, jj, kk, dtau, triad_mode, time_scheme)
end


function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat,
    triad_mode::Union{Triad2D, Triad3DIso},
    time_scheme::EulerMethod,
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int, nl_time_scale, dephasing_time, action_ref) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid

    (; action_rel_tol, increment_rel_tol, compute_dephasing_time) =
        state.namelists.triad

    # Reset the timescales so that an inactive physical cell does
    # not retain values calculated during an earlier timestep.
    nl_time_scale[ii, jj, kk] = Inf
    dephasing_time[ii, jj, kk] = Inf

    # Minimum significant spectral-cell action:
    #
    # A_floor = action_rel_tol * A_ref
    action_floor = action_rel_tol * action_ref[]

    #----------------------------------------------------------
    # Check whether the physical grid cell contains at least one
    # materially active spectral cell.
    #----------------------------------------------------------

    max_cell_action = 0.0

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        spectral_cell_width = delkp[kpi] * delm[mi]
        cell_action = wavespectrum[ii, jj, kk, kpi, mi] * spectral_cell_width

        if cell_action > max_cell_action
            max_cell_action = cell_action
        end
    end

    if max_cell_action <= action_floor
        return nothing
    end

    #----------------------------------------------------------
    # Compute the collision integral and nonlinear timescale
    # from the current wave spectrum.
    #----------------------------------------------------------

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    tau_nl = get_nl_time_scale(spec_tend, ii, jj, kk, action_rel_tol)
    nl_time_scale[ii, jj, kk] = tau_nl

    #----------------------------------------------------------
    # Maximum regularized relative action increment:
    #
    #                  Δτ |St| Δkp Δm
    # R = ------------------------------------------
    #     max(N Δkp Δm, A_floor)
    #----------------------------------------------------------

    max_relative_action_increment = 0.0

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]

        spectral_cell_width = delkp[kpi] * delm[mi]
        cell_action = was * spectral_cell_width
        predicted_action_increment = dtau * abs(st) * spectral_cell_width

        relative_action_increment =
            predicted_action_increment / max(cell_action, action_floor)

        if relative_action_increment > max_relative_action_increment
            max_relative_action_increment = relative_action_increment
        end
    end

    # The nonlinear interaction is insignificant throughout this
    # physical grid cell.
    if max_relative_action_increment <= increment_rel_tol
        return nothing
    end

    #----------------------------------------------------------
    # Compute the dephasing timescale only for a materially active
    # nonlinear interaction.
    #
    # The calculation is performed before the Euler update so that
    # wavespectrum, col_int, tau_nl and tau_pl all correspond to the
    # same pre-update state.
    #----------------------------------------------------------

    if compute_dephasing_time
       tau_pl = get_dephasing_time(state, ii, jj, kk, tau_nl, triad_mode)
        dephasing_time[ii, jj, kk] = tau_pl
    end

    #----------------------------------------------------------
    # Diagnostic: check whether the explicit Euler update would
    # produce negative wave-action density.
    #----------------------------------------------------------

    for mi in eachindex(m), kpi in eachindex(kp)
        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]

        was_new = was + dtau * st

        if was_new < 0.0
            relative_increment =
                was > 0.0 ? dtau * abs(st) / was : Inf

            println("")
            println("NEGATIVE WAVE SPECTRUM PREDICTED")
            println("  physical cell     = ", (ii, jj, kk))
            println("  spectral index    = ", (kpi, mi))
            println("  kp                = ", kp[kpi])
            println("  m                 = ", m[mi])
            println("  wavespectrum      = ", was)
            println("  col_int           = ", st)
            println("  dtau              = ", dtau)
            println("  dtau*|St|/N       = ", relative_increment)
            println("  tau_nl            = ", tau_nl)
            println("  dtau/tau_nl       = ",
                isfinite(tau_nl) ? dtau / tau_nl : 0.0)
            println("  predicted spectrum = ", was_new)
            println("")

            error("Negative wave-action density predicted by triad Euler update")
        end
    end

    #----------------------------------------------------------
    # Explicit Euler update.
    #----------------------------------------------------------

    @ivy for mi in eachindex(m), kpi in eachindex(kp)
        wavespectrum[ii, jj, kk, kpi, mi] +=
            dtau * col_int[ii, jj, kk, kpi, mi]
    end

    return nothing
end
#=
function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat, 
    triad_mode::Union{Triad2D, Triad3DIso},
    time_scheme::EulerMethod,
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int, nl_time_scale) = spec_tend
    (; kp, m) = spec_tend.spec_grid
    (; increment_rel_tol, action_abs_tol, action_rel_tol, st_abs_tol) = state.namelists.triad

    max_was = maximum(spec_tend.wavespectrum[ii, jj, kk, :, :])

    if max_was <= 1.0E-40
        #return if there is non significant wad in the physical grid cell
        return
    end

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)
    tau_nl = get_nl_time_scale(spec_tend, ii, jj, kk, action_abs_tol, action_rel_tol, st_abs_tol)
    nl_time_scale[ii, jj, kk] = tau_nl

    if (tau_nl * increment_rel_tol) > dtau
        #The nonlnear time scale is too large, interaction is not required in this grid cell
        return
    end 

    #Euler method 
    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)

        if  (dtau * abs(col_int[ii, jj, kk, kpi, mi])) < (increment_rel_tol * max_was) #cell wise exclusion for small col_int, the collision integral is too small just ignore it
            continue
        else
            wavespectrum[ii, jj, kk, kpi, mi] += dtau * col_int[ii, jj, kk, kpi, mi]  
        end
    end
    
end
=#
function update_wave_spectrum!(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    dtau::AbstractFloat, 
    triad_mode::Union{Triad2D, Triad3DIso},
    time_scheme::Rk2Step,
)
    (; spec_tend) = state.wkb
    (; wavespectrum, col_int) = spec_tend
    (; kp, m) = spec_tend.spec_grid

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    #RK2 method
    was_copy = spec_tend.wavespectrum
    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)

        if  col_int[kpi, mi] != 0
            
            wavespectrum[ii, jj, kk, kpi, mi] += 0.5 * dtau * col_int[kpi, mi]

        end

    end

    compute_scattering_integral!(state, ii, jj, kk, triad_mode)

    @ivy for mi in eachindex(m),
        kpi in eachindex(kp)
        
        if  col_int[kpi, mi] != 0
            wavespectrum[ii, jj, kk, kpi, mi] = was_copy[ii, jj, kk, kpi, mi] + dtau * col_int[kpi, mi]  
        end

    end

  
end
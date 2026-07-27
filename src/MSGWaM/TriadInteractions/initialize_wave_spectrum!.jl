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

function initialize_wave_spectrum!(state::State, 
    wkb_mode::Union{MultiColumn, SingleColumn, SteadyState}, 
    triad_mode::Union{Triad2D, Triad3DIso})

    (; master) = state.domain
    (; nthreads_triad, compute_dephasing_time, action_abs_tol, action_rel_tol, st_abs_tol) = state.namelists.triad
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; spec_tend) = state.wkb
    (; nl_time_scale) = spec_tend

    if master
        println(repeat("-", 80))
        println("")
        println("Triad interaction module activated with the model: ", triad_mode)
        println("Number of threads per rank: (nthreads_triad) = ", nthreads_triad)
        println("")
        println(repeat("-", 80))
    end

    get_wave_spectrum!(state, wkb_mode, triad_mode)

    nl_time_scale .= Inf

    @ivy for kk in k0:k1,
        jj in j0:j1,
        ii in i0:i1
        
            max_was = maximum(spec_tend.wavespectrum[ii, jj, kk, :, :])

            if max_was <= 1.0E-40
                #return if there is non significant wad in the physical grid cell
                return
            end

            compute_scattering_integral!(state, ii, jj, kk, triad_mode)
            tau_nl = get_nl_time_scale(spec_tend, ii, jj, kk, action_abs_tol, action_rel_tol, st_abs_tol)
            nl_time_scale[ii, jj, kk] = tau_nl        
    end

    if compute_dephasing_time
        spec_tend.consistency_time .= Inf
        compute_consistency_time!(state)
    end 

end


function apply_triad_interactions! end


function apply_triad_interactions!(state::State, dtau::AbstractFloat)
    (; wkb_mode) = state.namelists.wkb
    (; triad_mode) = state.namelists.triad
    apply_triad_interactions!(state, dtau, wkb_mode, triad_mode)
    return
end

function apply_triad_interactions!(state::State, 
    dtau::AbstractFloat, 
    wkb_mode::Union{NoWKB, MultiColumn, SingleColumn, SteadyState}, 
    triad_mode::NoTriad)
    return
end

function apply_triad_interactions!(state::State,
    dtau::AbstractFloat,
    wkb_mode::Union{MultiColumn, SingleColumn},
    triad_mode::Union{Triad2D, Triad3DIso})
    (; master) = state.domain

    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; spec_tend) = state.wkb
    (; nthreads_triad) = state.namelists.triad
    
    if Threads.nthreads() != nthreads_triad
        error("Julia started with Threads.nthreads()=$(Threads.nthreads()) but nthreads_triad=$nthreads_triad. Start Julia with --threads=$nthreads_triad (or set JULIA_NUM_THREADS).")
    end

    if master
        println(repeat("-", 80))
        println("Calling triad interaction module")
    end

    get_wave_spectrum!(state)
    wavespectrum_copy = deepcopy(spec_tend.wavespectrum)
    
    if master
        println("Updating wave action spectrum due to interactions")
    end

    @ivy for kk in k0:k1,
        jj in j0:j1,
        ii in i0:i1
        
            update_wave_spectrum!(state, ii, jj, kk, dtau, triad_mode)
        
    end
    
   
    get_ray_volumes!(state, wavespectrum_copy, triad_mode)

    if master
        println("Triad interaction module successfully called")
        println(repeat("-", 80))
    end

    
end




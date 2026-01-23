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

    (; domain, grid) = state
    (; branch) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = domain
    (; spec_tend) = state.wkb

    println("-------calling triad interaction module-------------")

    get_wave_spectrum!(state)
    wavespectrum_copy = deepcopy(spec_tend.wavespectrum)
   
    @ivy for kk in (k0-1):(k1+1),
        jj in (j0-1):(j1+1),
        ii in (i0-1):(i1+1)
    
        compute_scattering_integral!(state, ii, jj, kk, triad_mode)
        update_wave_spectrum!(state, ii, jj, kk, dtau, triad_mode)

    end
    
   
    get_ray_volumes!(state, wavespectrum_copy, triad_mode)

    println("----------Triad interaction module called successfully------------")

    
end




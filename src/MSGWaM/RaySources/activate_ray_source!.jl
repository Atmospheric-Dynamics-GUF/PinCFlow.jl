function activate_ray_source! end

function activate_ray_source!(state::State, rkstage::Integer)
    (; source_mode) = state.namelists.wkb

    activate_ray_source!(state, rkstage, source_mode)
    return
end

function activate_ray_source!(state::State, rkstage::Integer, source_mode::NoRaySource)
    return
end

function activate_ray_source!(state::State, rkstage::Integer, source_mode::OrographicSource)
    activate_orographic_source!(state)
    return
end

function activate_ray_source!(state::State, rkstage::Integer, source_mode::ContinuousSpectralSource)
    (; nstages) = state.time
    
    if rkstage == nstages
        activate_continuous_spectral_source!(state)
    end
    return
end
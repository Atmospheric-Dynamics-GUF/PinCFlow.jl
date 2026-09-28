function activate_ray_source! end

function activate_ray_source!(state::State)
    (; source_mode) = state.namelists.wkb

    activate_ray_source!(state, source_mode)
    return
end



function activate_ray_source!(state::State, source_mode::OrographicSource)
    activate_orographic_source!(state)
    return
end

function activate_ray_source!(state::State, source_mode::ContinuousSpectralSource)
    activate_continuous_spectral_source!(state)
    return
end
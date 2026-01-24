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
    get_wave_spectrum!(state, wkb_mode, triad_mode)
end


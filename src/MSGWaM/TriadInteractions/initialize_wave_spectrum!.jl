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
    (; nthreads_triad) = state.namelists.triad

    if master
        println(repeat("-", 80))
        println("")
        println("Triad interaction module activated with the model: ", triad_mode)
        println("Number of threads per rank: (nthreads_triad) = ", nthreads_triad)
        println("")
        println(repeat("-", 80))
    end

    get_wave_spectrum!(state, wkb_mode, triad_mode)
end


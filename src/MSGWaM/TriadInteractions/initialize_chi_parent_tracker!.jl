function initialize_chi_parent_tracker! end

function initialize_chi_parent_tracker!(state::State, mode_information)
    (; master) = state.domain
    (; wave_modes) = state.namelists.wkb
    (; chi_parent, chi_tracker) = state.spec_tend

    nidentified = length(mode_information)

    chi_parent .= NaN
    chi_tracker.ntracked[] = 0
    chi_tracker.kpi_ref .= 0
    chi_tracker.mi_ref .= 0
    chi_tracker.kpi_radius .= 0
    chi_tracker.mi_radius .= 0
    chi_tracker.valid .= false

    if nidentified > wave_modes
        master && @warn(
            "Initial spectral-mode identification found $nidentified components, " *
            "but wave_modes = $wave_modes. chi_parent tracking is disabled because " *
            "the initialized modes appear to have fragmented spectrally."
        )

        return nothing
    end

    chi_tracker.ntracked[] = nidentified

    for p in 1:nidentified
        info = mode_information[p]

        chi_tracker.kpi_ref[p] = info.peak_index[1]
        chi_tracker.mi_ref[p] = info.peak_index[2]

        nkp_mode = info.kpi_max - info.kpi_min + 1
        nm_mode = info.mi_max - info.mi_min + 1

        chi_tracker.kpi_radius[p] = max(1, cld(nkp_mode, 2) + 1)
        chi_tracker.mi_radius[p] = max(2, cld(nm_mode, 2) + 1)
    end

    if master && nidentified != wave_modes
        @warn(
            "$wave_modes wave modes were initialized, but only $nidentified " *
            "spectrally distinguishable components were identified. chi_parent " *
            "will track the identified components only; remaining entries stay NaN."
        )
    end

    return nothing
end
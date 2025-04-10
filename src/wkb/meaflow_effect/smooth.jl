function smooth!(state::State)
    (; lsmth_wkb) = state.namelists.wkb

    if lsmth_wkb
        smooth_tendencies!(state)
    end
    return
end

function smooth_tendencies!(state)
    (; sizex, sizey) = state.namelists.domain
    (; sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.integrals

    if sizey == 1
        if sizex > 1
            smooth_wkb!(dudt, state, sm_filter, XZ())
            smooth_wkb!(dvdt, state, sm_filter, XZ())
            smooth_wkb!(dthetadt, state, sm_filter, XZ())
        else
            error("Smoothing just in z not yet implemented.")
        end
    elseif sizex == 1
        smooth_wkb!(dudt, state, sm_filter, YZ())
        smooth_wkb!(dvdt, state, sm_filter, YZ())
        smooth_wkb!(dthetadt, state, sm_filter, YZ())
    elseif sizex > 1
        smooth_wkb!(dudt, state, sm_filter, XYZ())
        smooth_wkb!(dvdt, state, sm_filter, XYZ())
        smooth_wkb!(dthetadt, state, sm_filter, XYZ())
    end
    return
end

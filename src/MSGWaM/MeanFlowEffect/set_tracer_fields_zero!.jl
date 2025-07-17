function set_tracer_field_zero!(state)
    (; tracersetup) = state.namelists.tracer

    set_tracer_field_zero!(state, tracersetup)

    return
end

function set_tracer_field_zero!(state::State, tracersetup::NoTracer)
    return
end

function set_tracer_field_zero!(state::State, tracersetup::AbstractTracer)
    (; chiq0) = state.tracer.tracerforcings

    for field in fieldnames(TracerGWImpact)
        getfield(chiq0, field) .= 0.0 
    end
end
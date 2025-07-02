function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_vertical_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
            zboundaries,
            -,
        )
    end

    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
            zboundaries,
        )
    end

    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    tracersetup::AbstractTracer,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; tracerfluxes) = state.tracer

    if ko == 0
        for field in fieldnames(TracerFluxes)
            getfield(tracerfluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in fieldnames(TracerFluxes)
            getfield(tracerfluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

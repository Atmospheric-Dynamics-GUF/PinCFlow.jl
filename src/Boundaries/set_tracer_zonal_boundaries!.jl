function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerpredictands) = state.tracer

    for field in fieldnames(TracerPredictands)
        set_zonal_boundaries_of_field!(
            getfield(tracerpredictands, field),
            namelists,
            domain,
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; tracerreconstructions) = state.tracer

    for field in fieldnames(TracerReconstructions)
        set_zonal_boundaries_of_field!(
            getfield(tracerreconstructions, field),
            namelists,
            domain,
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_model::AbstractWKBMode,
    tracersetup::NoTracer,
)
    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:uchi, :vchi, :wchi)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::NoTracer,
)
    return 
end 

function set_tracer_zonal_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    wkb_mode::AbstractWKBMode,
    tracersetup::AbstractTracer,
)
    (; namelists, domain) = state
    (; chiq0) = state.tracer.tracerforcings

    for field in (:dchidt,)
        set_zonal_boundaries_of_field!(
            getfield(chiq0, field),
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    return
end
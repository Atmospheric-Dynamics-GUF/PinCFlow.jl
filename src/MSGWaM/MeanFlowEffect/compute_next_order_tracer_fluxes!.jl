function compute_next_order_tracer_fluxes! end

function compute_next_order_tracer_fluxes!(
    state::State,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    compute_next_order_tracer_fluxes!(
        state,
        tracer_setup,
        factor,
        r,
        i,
        j,
        k,
        xr,
        yr,
        zr,
        iray,
        jray,
        kray,
    )

    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    tracer_setup::NoTracer,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    return
end

function compute_next_order_tracer_fluxes!(
    state::State,
    tracer_setup::TracerOn,
    factor::AbstractFloat,
    r::Integer,
    i::Integer,
    j::Integer,
    k::Integer,
    xr::AbstractFloat,
    yr::AbstractFloat,
    zr::AbstractFloat,
    iray::Integer,
    jray::Integer,
    kray::Integer,
)
    return
end

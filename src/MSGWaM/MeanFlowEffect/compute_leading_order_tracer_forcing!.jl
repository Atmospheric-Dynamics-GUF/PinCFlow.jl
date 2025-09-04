function compute_leading_order_tracer_forcing!(
    state::State,
    ix::Integer,
    jy::Integer,
    kz::Integer,
    tracersetup::AbstractTracer,
)
    (; sizex, sizey) = state.namelists.domain
    (; dx, dy, dz) = state.grid
    (; uchi, vchi, wchi, dchidt) = state.tracer.tracerforcings.chiq0

    dchidt[ix, jy, kz] = 0.0

    if sizex > 1
        dchiu = (uchi[ix + 1, jy, kz] - uchi[ix - 1, jy, kz]) / (2.0 * dx)
    else
        dchiu = 0.0
    end

    if sizey > 1
        dchiv = (vchi[ix, jy + 1, kz] - vchi[ix, jy - 1, kz]) / (2.0 * dy)
    else
        dchiv = 0.0
    end

    dchiw = (wchi[ix, jy, kz + 1] - wchi[ix, jy, kz - 1]) / (2.0 * dz)

    dchidt[ix, jy, kz] = - (dchiu + dchiv + dchiw)

    return
end

function compute_leading_order_tracer_forcing!(
    state::State,
    ix::Integer,
    jy::Integer,
    kz::Integer,
    tracersetup::NoTracer,
)
    return
end

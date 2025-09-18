function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracersetup::AbstractTracer,
)
    (; sizex, sizey) = state.namelists.domain
    (; dx, dy, dz) = state.grid
    (; uchi, vchi, wchi, dchidt) = state.tracer.tracerforcings.chiq0

    @ivy dchidt[i, j, k] = 0.0

    @ivy if sizex > 1
        dchiu = (uchi[i + 1, j, k] - uchi[i - 1, j, k]) / (2.0 * dx)
    else
        dchiu = 0.0
    end

    @ivy if sizey > 1
        dchiv = (vchi[i, j + 1, k] - vchi[i, j - 1, k]) / (2.0 * dy)
    else
        dchiv = 0.0
    end

    @ivy dchiw = (wchi[i, j, k + 1] - wchi[i, j, k - 1]) / (2.0 * dz)

    @ivy dchidt[i, j, k] = -(dchiu + dchiv + dchiw)

    return
end

function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracersetup::NoTracer,
)
    return
end

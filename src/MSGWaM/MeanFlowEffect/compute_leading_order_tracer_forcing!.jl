"""
```julia
compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::AbstractTracer,
)
```

Compute and return the leading-order tracer forcing at ``\\left(i, j, k\\right)``.

```julia
compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `tracer_setup`: General tracer-transport configuration.
"""
function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::AbstractTracer,
)
    (; x_size, y_size) = state.namelists.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; uchi, vchi, wchi, dchidt) = state.tracer.tracerforcings.chiq0

    @ivy dchidt[i, j, k] = 0.0

    @ivy if x_size > 1
        dchiu =
            (uchi[i + 1, j, k] - uchi[i - 1, j, k]) / (2.0 * dx) +
            met[i, j, k, 1, 3] * (uchi[i, j, k + 1] - uchi[i, j, k - 1]) /
            (2.0 * dz)
    else
        dchiu = 0.0
    end

    @ivy if y_size > 1
        dchiv =
            (vchi[i, j + 1, k] - vchi[i, j - 1, k]) / (2.0 * dy) +
            met[i, j, k, 2, 3] * (vchi[i, j, k + 1] - vchi[i, j, k - 1]) /
            (2.0 * dz)
    else
        dchiv = 0.0
    end

    @ivy dchiw =
        (wchi[i, j, k + 1] - wchi[i, j, k - 1]) / (2.0 * jac[i, j, k] * dz)

    @ivy dchidt[i, j, k] = -(dchiu + dchiv + dchiw)

    return
end

function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::NoTracer,
)
    return
end

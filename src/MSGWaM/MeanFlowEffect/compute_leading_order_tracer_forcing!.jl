"""
```julia
compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
```

Compute the leading-order tracer forcing by dispatching to the tracer-setup-specific method.

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

```julia
compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::TracerOn,
)
```

Compute the leading-order tracer forcing at ``\\left(i, j, k\\right)``.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `tracer_setup`: General tracer-transport configuration.
"""
function compute_leading_order_tracer_forcing! end

function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    compute_leading_order_tracer_forcing!(state, i, j, k, tracer_setup)
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

function compute_leading_order_tracer_forcing!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::TracerOn,
)
    (; x_size, y_size) = state.namelists.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; uchi0, vchi0, wchi0) = state.tracer.tracerwkbintegrals
    (; dchidt0) = state.tracer.tracerwkbtendencies
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    @ivy dchidt0[i, j, k] = 0.0

    @ivy if x_size > 1
        dchiu =
            (uchi0[i + 1, j, k] - uchi0[i - 1, j, k]) / (2.0 * dx) +
            met[i, j, k, 1, 3] * (uchi0[i, j, k + 1] - uchi0[i, j, k - 1]) /
            (2.0 * dz)
    else
        dchiu = 0.0
    end

    @ivy if y_size > 1
        dchiv =
            (vchi0[i, j + 1, k] - vchi0[i, j - 1, k]) / (2.0 * dy) +
            met[i, j, k, 2, 3] * (vchi0[i, j, k + 1] - vchi0[i, j, k - 1]) /
            (2.0 * dz)
    else
        dchiv = 0.0
    end

    @ivy dchiw =
        (wchi0[i, j, k + 1] - wchi0[i, j, k - 1]) / (2.0 * jac[i, j, k] * dz)

    @ivy dchidt0[i, j, k] =
        -(rho[i, j, k] + rhobar[i, j, k]) / rhobar[i, j, k] *
        (dchiu + dchiv + dchiw)

    return
end

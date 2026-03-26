"""
```julia
compute_gw_tracer_tendencies!(state::State, i::Integer, j::Integer, k::Integer)
```

Compute the leading-order tracer forcing by dispatching to the tracer-setup-specific method.

```julia
compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_tendencies!(
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
function compute_gw_tracer_tendencies! end

function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; tracer_setup) = state.namelists.tracer

    compute_gw_tracer_tendencies!(state, i, j, k, tracer_setup)
    return
end

function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::NoTracer,
)
    return
end

function compute_gw_tracer_tendencies!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    tracer_setup::TracerOn,
)
    (; x_size, y_size) = state.namelists.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; uchi0, vchi0, wchi0, uchi1, vchi1, wchi1, qchi) =
        state.tracer.tracerwkbintegrals
    (; dchidt0, dchidtq, dchidt1) = state.tracer.tracerwkbtendencies
    (; rho) = state.variables.predictands
    (; rhobar, thetabar) = state.atmosphere

    @ivy dchidt0[i, j, k] = 0.0
    @ivy dchidtq[i, j, k] = 0.0
    @ivy dchidt1[i, j, k] = 0.0

    @ivy if x_size > 1
        dchiu =
            (uchi0[i + 1, j, k] - uchi0[i - 1, j, k]) / (2.0 * dx) +
            met[i, j, k, 1, 3] * (uchi0[i, j, k + 1] - uchi0[i, j, k - 1]) /
            (2.0 * dz)
        dchiu1 =
            (
                rhobar[i + 1, j, k] *
                thetabar[i + 1, j, k] *
                uchi1[i + 1, j, k] -
                rhobar[i - 1, j, k] *
                thetabar[i - 1, j, k] *
                uchi1[i - 1, j, k]
            ) / (2.0 * dx) +
            met[i, j, k, 1, 3] * (
                rhobar[i, j, k + 1] *
                thetabar[i, j, k + 1] *
                uchi1[i, j, k + 1] -
                rhobar[i, j, k - 1] *
                thetabar[i, j, k - 1] *
                uchi1[i, j, k - 1]
            ) / (2.0 * dz)
    else
        dchiu = 0.0
        dchiu1 = 0.0
    end

    @ivy if y_size > 1
        dchiv =
            (vchi0[i, j + 1, k] - vchi0[i, j - 1, k]) / (2.0 * dy) +
            met[i, j, k, 2, 3] * (vchi0[i, j, k + 1] - vchi0[i, j, k - 1]) /
            (2.0 * dz)
        dchiv1 =
            (
                rhobar[i, j + 1, k] *
                thetabar[i, j + 1, k] *
                vchi1[i, j + 1, k] -
                rhobar[i, j - 1, k] *
                thetabar[i, j - 1, k] *
                vchi1[i, j - 1, k]
            ) / (2.0 * dy) +
            met[i, j, k, 2, 3] * (
                rhobar[i, j, k + 1] *
                thetabar[i, j, k + 1] *
                vchi1[i, j, k + 1] -
                rhobar[i, j, k - 1] *
                thetabar[i, j, k - 1] *
                vchi1[i, j, k - 1]
            ) / (2.0 * dz)
    else
        dchiv = 0.0
        dchiv1 = 0.0
    end

    @ivy dchiw =
        (wchi0[i, j, k + 1] - wchi0[i, j, k - 1]) / (2.0 * jac[i, j, k] * dz)
    @ivy dchiw1 =
        (
            rhobar[i, j, k + 1] * thetabar[i, j, k + 1] * wchi1[i, j, k + 1] -
            rhobar[i, j, k - 1] * thetabar[i, j, k - 1] * wchi1[i, j, k - 1]
        ) / (2.0 * jac[i, j, k] * dz)

    @ivy dchidt0[i, j, k] =
        -(rho[i, j, k] + rhobar[i, j, k]) / rhobar[i, j, k] *
        (dchiu + dchiv + dchiw)

    @ivy dchidt1[i, j, k] =
        -(rho[i, j, k] + rhobar[i, j, k]) /
        (2 * rhobar[i, j, k] * thetabar[i, j, k]) * (dchiu1 + dchiv1 + dchiw1)

    @ivy dchidtq[i, j, k] =
        -(rho[i, j, k] + rhobar[i, j, k]) / 2 *
        (qchi[i, j, k + 1] - qchi[i, j, k - 1]) / (2.0 * jac[i, j, k] * dz)

    return
end

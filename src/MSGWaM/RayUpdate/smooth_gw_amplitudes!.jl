"""
```julia
smooth_gw_amplitudes!(state::State)
```

Apply spatial smoothing to gravity-wave tendency fields by dispatching to a method specific for the chosen filter (`state.namelists.wkb.filter_type`) and dimensionality of the domain.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::BoxFilter,
    direction::XYZ,
)
```

Apply a 3D box filter to smooth in all spatial directions.

Applies the moving average

```math
\\tilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 3} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::BoxFilter,
    direction::XZ,
)
```

Apply a 2D box filter to smooth in ``\\hat{x}`` and ``\\hat{z}``.

Applies the moving average

```math
\\tilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::BoxFilter,
    direction::YZ,
)
```

Apply a 2D box filter to smooth in ``\\hat{y}`` and ``\\hat{z}``.

Applies the moving average

```math
\\tilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::BoxFilter,
    direction::Z,
)
```

Apply a 1D box filter to smooth in ``\\hat{z}``.

Applies the moving average

```math
\\tilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 1} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::XYZ,
)
```

Apply a 3D Shapiro filter to smooth in all spatial directions.

A 1D Shapiro filter is applied sequentially in ``\\hat{x}``, ``\\hat{y}`` and ``\\hat{z}``.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::XZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\hat{x}`` and ``\\hat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\hat{x}`` and ``\\hat{z}``.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::YZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\hat{y}`` and ``\\hat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\hat{y}`` and ``\\hat{z}``.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::Z,
)
```

Apply a 1D Shapiro filter to smooth in ``\\hat{z}``.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::Y,
)
```

Apply a 1D Shapiro filter to smooth in ``\\hat{y}``.

```julia
smooth_gw_amplitudes!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::ShapiroFilter,
    direction::X,
)
```

Apply a 1D Shapiro filter to smooth in ``\\hat{x}``.

```julia
smooth_gw_amplitudes!(state::State, tracer_setup::TracerOn)
```

Apply smoothing to tracer tendencies.

```julia
smooth_gw_amplitudes!(state::State, tracer_setup::NoTracer)
```

Return for configurations without tracer transport.

# Arguments

  - `state`: Model state.

  - `output`: Field to smooth.

  - `filter_type`: Filter type.

  - `direction`: Directions to smooth in.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.MSGWaM.MeanFlowEffect.apply_shapiro_filter!`](@ref)
"""
function smooth_gw_amplitudes! end

function smooth_gw_amplitudes!(state::State)
    (; tracer_setup) = state.namelists.tracer
    (; smooth_tendencies) = state.namelists.wkb

    if !smooth_tendencies
        return
    end

    smooth_gw_amplitudes!(state, tracer_setup)

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::BoxFilter,
    direction::XYZ,
) where {T <: Real}
    (; nbx, nby, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_amplitudes!: nbx < filter_order!")
    end
    if nby < filter_order
        error("Error in smooth_gw_amplitudes!: nby < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_amplitudes!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[
                    (i - filter_order):(i + filter_order),
                    (j - filter_order):(j + filter_order),
                    (k - filter_order):(k + filter_order),
                ],
            ) / (2 * filter_order + 1)^3
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::BoxFilter,
    direction::XZ,
) where {T <: Real}
    (; nbx, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_amplitudes!: nbx < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_amplitudes!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[
                    (i - filter_order):(i + filter_order),
                    j,
                    (k - filter_order):(k + filter_order),
                ],
            ) / (2 * filter_order + 1)^2
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::BoxFilter,
    direction::YZ,
) where {T <: Real}
    (; nby, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < filter_order
        error("Error in smooth_gw_amplitudes!: nby < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_amplitudes!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[
                    i,
                    (j - filter_order):(j + filter_order),
                    (k - filter_order):(k + filter_order),
                ],
            ) / (2 * filter_order + 1)^2
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::BoxFilter,
    direction::Z,
) where {T <: Real}
    (; nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < filter_order
        error("Error in smooth_gw_amplitudes!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(input[i, j, (k - filter_order):(k + filter_order)]) /
            (2 * filter_order + 1)
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::XYZ,
) where {T <: Real}
    smooth_gw_amplitudes!(output, state, filter_type, X())
    smooth_gw_amplitudes!(output, state, filter_type, Y())
    smooth_gw_amplitudes!(output, state, filter_type, Z())
    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::XZ,
) where {T <: Real}
    smooth_gw_amplitudes!(output, state, filter_type, X())
    smooth_gw_amplitudes!(output, state, filter_type, Z())
    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::YZ,
) where {T <: Real}
    smooth_gw_amplitudes!(output, state, filter_type, Y())
    smooth_gw_amplitudes!(output, state, filter_type, Z())
    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::Z,
) where {T <: Real}
    (; nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nxx, nyy, k0, k1) = state.domain

    if nbz < filter_order
        error("Error in smooth_gw_amplitudes!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for j in 1:nyy, i in 1:nxx
        apply_shapiro_filter!(
            output[i, j, :],
            input[i, j, :],
            k0:k1,
            Val(filter_order),
        )
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::Y,
) where {T <: Real}
    (; nby) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nxx, nzz, j0, j1) = state.domain

    if nby < filter_order
        error("Error in smooth_gw_amplitudes!: nby < filter_order!")
    end

    input = copy(output)
    @ivy for k in 1:nzz, i in 1:nxx
        apply_shapiro_filter!(
            output[i, :, k],
            input[i, :, k],
            j0:j1,
            Val(filter_order),
        )
    end

    return
end

function smooth_gw_amplitudes!(
    output::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    state::State,
    filter_type::ShapiroFilter,
    direction::X,
) where {T <: Real}
    (; nbx) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nyy, nzz, i0, i1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_amplitudes!: nbx < filter_order!")
    end

    input = copy(output)
    @ivy for k in 1:nzz, j in 1:nyy
        apply_shapiro_filter!(
            output[:, j, k],
            input[:, j, k],
            i0:i1,
            Val(filter_order),
        )
    end

    return
end

function smooth_gw_amplitudes!(state::State, tracer_setup::TracerOn)
    (; x_size, y_size) = state.namelists.domain
    (; smooth_tendencies, filter_type) = state.namelists.wkb
    (; uhat, vhat, what, bhat, pihat, chihat) = state.wkb.integrals

    if !smooth_tendencies
        return
    end

    if x_size == y_size == 1
        smooth_gw_amplitudes!(uhat, state, filter_type, Z())
        smooth_gw_amplitudes!(vhat, state, filter_type, Z())
        smooth_gw_amplitudes!(what, state, filter_type, Z())
        smooth_gw_amplitudes!(bhat, state, filter_type, Z())
        smooth_gw_amplitudes!(pihat, state, filter_type, Z())
        smooth_gw_amplitudes!(chihat, state, filter_type, Z())
    elseif x_size == 1
        smooth_gw_amplitudes!(uhat, state, filter_type, YZ())
        smooth_gw_amplitudes!(vhat, state, filter_type, YZ())
        smooth_gw_amplitudes!(what, state, filter_type, YZ())
        smooth_gw_amplitudes!(bhat, state, filter_type, YZ())
        smooth_gw_amplitudes!(pihat, state, filter_type, YZ())
        smooth_gw_amplitudes!(chihat, state, filter_type, YZ())
    elseif y_size == 1
        smooth_gw_amplitudes!(uhat, state, filter_type, XZ())
        smooth_gw_amplitudes!(vhat, state, filter_type, XZ())
        smooth_gw_amplitudes!(what, state, filter_type, XZ())
        smooth_gw_amplitudes!(bhat, state, filter_type, XZ())
        smooth_gw_amplitudes!(pihat, state, filter_type, XZ())
        smooth_gw_amplitudes!(chihat, state, filter_type, XZ())
    else
        smooth_gw_amplitudes!(uhat, state, filter_type, XYZ())
        smooth_gw_amplitudes!(vhat, state, filter_type, XYZ())
        smooth_gw_amplitudes!(what, state, filter_type, XYZ())
        smooth_gw_amplitudes!(bhat, state, filter_type, XYZ())
        smooth_gw_amplitudes!(pihat, state, filter_type, XYZ())
        smooth_gw_amplitudes!(chihat, state, filter_type, XYZ())
    end

    return
end

function smooth_gw_amplitudes!(state::State, tracer_setup::NoTracer)
    return
end

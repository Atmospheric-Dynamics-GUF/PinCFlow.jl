"""
```julia
smooth_gw_tendencies!(state::State)
```

Apply spatial smoothing to gravity-wave tendency fields by dispatching to a method specific for the chosen filter (`state.namelists.wkb.filter_type`) and dimensionality of the domain.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::XYZ,
)
```

Apply a 3D box filter to smooth in all spatial directions.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 3} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::XZ,
)
```

Apply a 2D box filter to smooth in ``\\widehat{x}`` and ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::YZ,
)
```

Apply a 2D box filter to smooth in ``\\widehat{y}`` and ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::Z,
)
```

Apply a 1D box filter to smooth in ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 1} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.filter_order`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::XYZ,
)
```

Apply a 3D Shapiro filter to smooth in all spatial directions.

A 1D Shapiro filter is applied sequentially in ``\\widehat{x}``, ``\\widehat{y}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::XZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\widehat{x}`` and ``\\widehat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\widehat{x}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::YZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\widehat{y}`` and ``\\widehat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\widehat{y}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::Z,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::Y,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{y}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::X,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{x}``.

```julia
smooth_gw_tendencies!(state::State, tracer_setup::TracerOn)
```

Apply smoothing to tracer tendencies.

```julia
smooth_gw_tendencies!(state::State, tracer_setup::NoTracer)
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
function smooth_gw_tendencies! end

function smooth_gw_tendencies!(state::State)
    (; x_size, y_size) = state.namelists.domain
    (; smooth_tendencies, filter_type) = state.namelists.wkb
    (; dudt, dvdt, dthetadt, shear) = state.wkb.tendencies
    (; tracer_setup) = state.namelists.tracer

    if !smooth_tendencies
        return
    end

    if x_size == y_size == 1
        smooth_gw_tendencies!(dudt, state, filter_type, Z())
        smooth_gw_tendencies!(dvdt, state, filter_type, Z())
        smooth_gw_tendencies!(dthetadt, state, filter_type, Z())
        smooth_gw_tendencies!(shear, state, filter_type, Z())
    elseif x_size == 1
        smooth_gw_tendencies!(dudt, state, filter_type, YZ())
        smooth_gw_tendencies!(dvdt, state, filter_type, YZ())
        smooth_gw_tendencies!(dthetadt, state, filter_type, YZ())
        smooth_gw_tendencies!(shear, state, filter_type, YZ())
    elseif y_size == 1
        smooth_gw_tendencies!(dudt, state, filter_type, XZ())
        smooth_gw_tendencies!(dvdt, state, filter_type, XZ())
        smooth_gw_tendencies!(dthetadt, state, filter_type, XZ())
        smooth_gw_tendencies!(shear, state, filter_type, XZ())
    else
        smooth_gw_tendencies!(dudt, state, filter_type, XYZ())
        smooth_gw_tendencies!(dvdt, state, filter_type, XYZ())
        smooth_gw_tendencies!(dthetadt, state, filter_type, XYZ())
        smooth_gw_tendencies!(shear, state, filter_type, XYZ())
    end

    smooth_gw_tendencies!(state, tracer_setup)

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::XYZ,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_tendencies!: nbx < filter_order!")
    end
    if nby < filter_order
        error("Error in smooth_gw_tendencies!: nby < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_tendencies!: nbz < filter_order!")
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

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::XZ,
)
    (; nbx, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_tendencies!: nbx < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_tendencies!: nbz < filter_order!")
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

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::YZ,
)
    (; nby, nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < filter_order
        error("Error in smooth_gw_tendencies!: nby < filter_order!")
    end
    if nbz < filter_order
        error("Error in smooth_gw_tendencies!: nbz < filter_order!")
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

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Box,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < filter_order
        error("Error in smooth_gw_tendencies!: nbz < filter_order!")
    end

    input = copy(output)
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(input[i, j, (k - filter_order):(k + filter_order)]) /
            (2 * filter_order + 1)
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::XYZ,
)
    smooth_gw_tendencies!(output, state, filter_type, X())
    smooth_gw_tendencies!(output, state, filter_type, Y())
    smooth_gw_tendencies!(output, state, filter_type, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::XZ,
)
    smooth_gw_tendencies!(output, state, filter_type, X())
    smooth_gw_tendencies!(output, state, filter_type, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::YZ,
)
    smooth_gw_tendencies!(output, state, filter_type, Y())
    smooth_gw_tendencies!(output, state, filter_type, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nxx, nyy, k0, k1) = state.domain

    if nbz < filter_order
        error("Error in smooth_gw_tendencies!: nbz < filter_order!")
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

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::Y,
)
    (; nby) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nxx, nzz, j0, j1) = state.domain

    if nby < filter_order
        error("Error in smooth_gw_tendencies!: nby < filter_order!")
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

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    filter_type::Shapiro,
    direction::X,
)
    (; nbx) = state.namelists.domain
    (; filter_order) = state.namelists.wkb
    (; nyy, nzz, i0, i1) = state.domain

    if nbx < filter_order
        error("Error in smooth_gw_tendencies!: nbx < filter_order!")
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

function smooth_gw_tendencies!(state::State, tracer_setup::TracerOn)
    (; x_size, y_size) = state.namelists.domain
    (; smooth_tendencies, filter_type) = state.namelists.wkb
    (; dchidt) = state.tracer.tracerforcings.chiq0

    if !smooth_tendencies
        return
    end

    if x_size == y_size == 1
        smooth_gw_tendencies!(dchidt, state, filter_type, Z())
    elseif x_size == 1
        smooth_gw_tendencies!(dchidt, state, filter_type, YZ())
    elseif y_size == 1
        smooth_gw_tendencies!(dchidt, state, filter_type, XZ())
    else
        smooth_gw_tendencies!(dchidt, state, filter_type, XYZ())
    end

    return
end

function smooth_gw_tendencies!(state::State, tracer_setup::NoTracer)
    return
end

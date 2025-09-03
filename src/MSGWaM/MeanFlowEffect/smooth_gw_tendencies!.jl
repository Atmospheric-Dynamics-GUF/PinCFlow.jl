"""
```julia
smooth_gw_tendencies!(state::State)
```

Apply spatial smoothing to gravity-wave tendency fields by dispatching to a method specific for the chosen filter (`state.namelists.wkb.sm_filter`) and dimensionality of the domain.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XYZ,
)
```

Apply a 3D box filter to smooth in all spatial directions.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 3} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.nsmth_wkb`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XZ,
)
```

Apply a 2D box filter to smooth in ``\\widehat{x}`` and ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\lambda = i - N_\\mathrm{s}}^{i + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{\\lambda, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.nsmth_wkb`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::YZ,
)
```

Apply a 2D box filter to smooth in ``\\widehat{y}`` and ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 2} \\sum\\limits_{\\mu = j - N_\\mathrm{s}}^{j + N_\\mathrm{s}} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, \\mu, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.nsmth_wkb`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::Z,
)
```

Apply a 1D box filter to smooth in ``\\widehat{z}``.

Applies the moving average

```math
\\widetilde{\\phi}_{i, j, k} = \\left(2 N_\\mathrm{s} + 1\\right)^{- 1} \\sum\\limits_{\\nu = k - N_\\mathrm{s}}^{k + N_\\mathrm{s}} \\phi_{i, j, \\nu},
```

where ``N_\\mathrm{s}`` is the order of the filter (`state.namelists.wkb.nsmth_wkb`).

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XYZ,
)
```

Apply a 3D Shapiro filter to smooth in all spatial directions.

A 1D Shapiro filter is applied sequentially in ``\\widehat{x}``, ``\\widehat{y}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\widehat{x}`` and ``\\widehat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\widehat{x}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::YZ,
)
```

Apply a 2D Shapiro filter to smooth in ``\\widehat{y}`` and ``\\widehat{z}``.

A 1D Shapiro filter is applied sequentially in ``\\widehat{y}`` and ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Z,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{z}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Y,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{y}``.

```julia
smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::X,
)
```

Apply a 1D Shapiro filter to smooth in ``\\widehat{x}``.

# Arguments

  - `state`: Model state.

  - `output`: Field to smooth.

  - `sm_filter`: Filter type.

  - `direction`: Directions to smooth in.

# See also

  - [`PinCFlow.MSGWaM.MeanFlowEffect.apply_shapiro_filter!`](@ref)
"""
function smooth_gw_tendencies! end

function smooth_gw_tendencies!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies
    (; tracersetup) = state.namelists.tracer

    if !lsmth_wkb
        return
    end

    if sizex == sizey == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, Z())
        smooth_gw_tendencies!(dvdt, state, sm_filter, Z())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, Z())
    elseif sizex == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, YZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, YZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, YZ())
    elseif sizey == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, XZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, XZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, XZ())
    else
        smooth_gw_tendencies!(dudt, state, sm_filter, XYZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, XYZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, XYZ())
    end

    smooth_gw_tendencies!(state, tracersetup)

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XYZ,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end
    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    (i - nsmth_wkb):(i + nsmth_wkb),
                    (j - nsmth_wkb):(j + nsmth_wkb),
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^3
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XZ,
)
    (; nbx, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    (i - nsmth_wkb):(i + nsmth_wkb),
                    j,
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^2
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::YZ,
)
    (; nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    i,
                    (j - nsmth_wkb):(j + nsmth_wkb),
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^2
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(input[i, j, (k - nsmth_wkb):(k + nsmth_wkb)]) /
            (2 * nsmth_wkb + 1)
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XYZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, X())
    smooth_gw_tendencies!(output, state, sm_filter, Y())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, X())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::YZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, Y())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nxx, nyy, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for j in 1:nyy, i in 1:nxx
        @views apply_shapiro_filter!(
            output[i, j, :],
            input[i, j, :],
            (k0, k1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Y,
)
    (; nby) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nxx, nzz, j0, j1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end

    input = copy(output)
    for k in 1:nzz, i in 1:nxx
        @views apply_shapiro_filter!(
            output[i, :, k],
            input[i, :, k],
            (j0, j1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::X,
)
    (; nbx) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nyy, nzz, i0, i1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end

    input = copy(output)
    for k in 1:nzz, j in 1:nyy
        @views apply_shapiro_filter!(
            output[:, j, k],
            input[:, j, k],
            (i0, i1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_gw_tendencies!(state::State, tracersetup::AbstractTracer)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dchidt) = state.tracer.tracerforcings.chiq0

    if !lsmth_wkb
        return
    end

    if sizex == sizey == 1
        smooth_gw_tendencies!(dchidt, state, sm_filter, Z())
    elseif sizex == 1
        smooth_gw_tendencies!(dchidt, state, sm_filter, YZ())
    elseif sizey == 1
        smooth_gw_tendencies!(dchidt, state, sm_filter, XZ())
    else
        smooth_gw_tendencies!(dchidt, state, sm_filter, XYZ())
    end

    return
end

function smooth_gw_tendencies!(state::State, tracersetup::NoTracer)
    return
end
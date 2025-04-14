function smooth_tendencies!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.integrals

    if !lsmth_wkb
        return
    end

    if sizex == sizey == 1
        smooth_tendencies!(dudt, state, sm_filter, Z())
        smooth_tendencies!(dvdt, state, sm_filter, Z())
        smooth_tendencies!(dthetadt, state, sm_filter, Z())
    elseif sizex == 1
        smooth_tendencies!(dudt, state, sm_filter, YZ())
        smooth_tendencies!(dvdt, state, sm_filter, YZ())
        smooth_tendencies!(dthetadt, state, sm_filter, YZ())
    elseif sizey == 1
        smooth_tendencies!(dudt, state, sm_filter, XZ())
        smooth_tendencies!(dvdt, state, sm_filter, XZ())
        smooth_tendencies!(dthetadt, state, sm_filter, XZ())
    else
        smooth_tendencies!(dudt, state, sm_filter, XYZ())
        smooth_tendencies!(dvdt, state, sm_filter, XYZ())
        smooth_tendencies!(dthetadt, state, sm_filter, XYZ())
    end

    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XYZ,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_tendencies!: nbx < nsmth_wkb!")
    end
    if nby < nsmth_wkb
        error("Error in smooth_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_tendencies!: nbz < nsmth_wkb!")
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

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XZ,
)
    (; nbx, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_tendencies!: nbx < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_tendencies!: nbz < nsmth_wkb!")
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

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::YZ,
)
    (; nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_tendencies!: nbz < nsmth_wkb!")
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

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(input[i, j, (k - nsmth_wkb):(k + nsmth_wkb)]) /
            (2 * nsmth_wkb + 1)
    end

    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XYZ,
)
    smooth_tendencies!(output, state, sm_filter, X())
    smooth_tendencies!(output, state, sm_filter, Y())
    smooth_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XZ,
)
    smooth_tendencies!(output, state, sm_filter, X())
    smooth_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::YZ,
)
    smooth_tendencies!(output, state, sm_filter, Y())
    smooth_tendencies!(output, state, sm_filter, Z())
    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for j in j0:j1, i in i0:i1
        @views smooth_tendencies!(
            input[i, j, :],
            output[i, j, :],
            (k0, k1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Y,
)
    (; nby) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_tendencies!: nby < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, i in i0:i1
        @views smooth_tendencies!(
            input[i, :, k],
            output[i, :, k],
            (j0, j1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::X,
)
    (; nbx) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_tendencies!: nbx < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1
        @views smooth_tendencies!(
            input[:, j, k],
            output[:, j, k],
            (i0, i1),
            Val(nsmth_wkb),
        )
    end

    return
end

function smooth_tendencies!(
    input::AbstractVector{<:AbstractFloat},
    output::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{1},
)
    for i in bounds[1]:bounds[2]
        output[i] = (input[i - 1] + input[i + 1] + 2 * input[i]) / 4
    end
    return
end

function smooth_tendencies!(
    input::AbstractVector{<:AbstractFloat},
    output::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{2},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                -input[i - 2] - input[i + 2] +
                4 * (input[i - 1] + input[i + 1]) +
                10 * input[i]
            ) / 16
    end
    return
end

function smooth_tendencies!(
    input::AbstractVector{<:AbstractFloat},
    output::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{3},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                input[i - 3] + input[i + 3] -
                6 * (input[i - 2] + input[i + 2]) +
                15 * (input[i - 1] + input[i + 1]) +
                44 * input[i]
            ) / 64
    end
    return
end

function smooth_tendencies!(
    input::AbstractVector{<:AbstractFloat},
    output::AbstractVector{<:AbstractFloat},
    bounds::NTuple{2, <:Integer},
    order::Val{4},
)
    for i in bounds[1]:bounds[2]
        output[i] =
            (
                -input[i - 4] - input[i + 4] +
                8 * (input[i - 3] + input[i + 3]) -
                28 * (input[i - 2] + input[i + 2]) +
                56 * (input[i - 1] + input[i + 1]) +
                186 * input[i]
            ) / 256
    end
    return
end
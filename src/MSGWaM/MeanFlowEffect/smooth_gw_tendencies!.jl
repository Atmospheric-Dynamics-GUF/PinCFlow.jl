function smooth_gw_tendencies!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

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

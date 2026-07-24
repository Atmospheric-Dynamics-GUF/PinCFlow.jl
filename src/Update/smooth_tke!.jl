function smooth_tke! end

function smooth_tke!(state::State)
    (; x_size, y_size) = state.namelists.domain
    (; smooth_turbulence, turbulence_filter_type) = state.namelists.turbulence
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere

    if !smooth_turbulence
        return
    end

    tke ./= (rho .+ rhobar)

    @dispatch_turbulence_filter_type if x_size == y_size == 1
        smooth_tke!(tke, state, Val(turbulence_filter_type), Z())
    elseif x_size == 1
        smooth_tke!(tke, state, Val(turbulence_filter_type), YZ())
    elseif y_size == 1
        smooth_tke!(tke, state, Val(turbulence_filter_type), XZ())
    else
        smooth_tke!(tke, state, Val(turbulence_filter_type), XYZ())
    end

    tke .*= (rho .+ rhobar)

    return
end

@ivy function smooth_tke!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    turbulence_filter_type::Val{:BoxFilter},
    direction::XYZ,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; turbulence_filter_order) = state.namelists.turbulence
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid

    if nbx < turbulence_filter_order
        error("Error in smooth_tke!: nbx < turbulence_filter_order!")
    end
    if nby < turbulence_filter_order
        error("Error in smooth_tke!: nby < turbulence_filter_order!")
    end
    if nbz < turbulence_filter_order
        error("Error in smooth_tke!: nbz < turbulence_filter_order!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[ii, jj, kk] * jac[ii, jj, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                jj in
                (j - turbulence_filter_order):(j + turbulence_filter_order),
                ii in
                (i - turbulence_filter_order):(i + turbulence_filter_order)
            ) / sum(
                jac[ii, jj, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                jj in
                (j - turbulence_filter_order):(j + turbulence_filter_order),
                ii in
                (i - turbulence_filter_order):(i + turbulence_filter_order)
            )
    end

    return
end

@ivy function smooth_tke!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    turbulence_filter_type::Val{:BoxFilter},
    direction::XZ,
)
    (; nbx, nbz) = state.namelists.domain
    (; turbulence_filter_order) = state.namelists.turbulence
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid

    if nbx < turbulence_filter_order
        error("Error in smooth_tke!: nbx < turbulence_filter_order!")
    end
    if nbz < turbulence_filter_order
        error("Error in smooth_tke!: nbz < turbulence_filter_order!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[ii, j, kk] * jac[ii, j, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                ii in
                (i - turbulence_filter_order):(i + turbulence_filter_order)
            ) / sum(
                jac[ii, j, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                ii in
                (i - turbulence_filter_order):(i + turbulence_filter_order)
            )
    end

    return
end

@ivy function smooth_tke!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    turbulence_filter_type::Val{:BoxFilter},
    direction::YZ,
)
    (; nby, nbz) = state.namelists.domain
    (; turbulence_filter_order) = state.namelists.turbulence
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid

    if nby < turbulence_filter_order
        error("Error in smooth_tke!: nby < turbulence_filter_order!")
    end
    if nbz < turbulence_filter_order
        error("Error in smooth_tke!: nbz < turbulence_filter_order!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[i, jj, kk] * jac[i, jj, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                jj in
                (j - turbulence_filter_order):(j + turbulence_filter_order)
            ) / sum(
                jac[i, jj, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order),
                jj in
                (j - turbulence_filter_order):(j + turbulence_filter_order)
            )
    end

    return
end

@ivy function smooth_tke!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    turbulence_filter_type::Val{:BoxFilter},
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; turbulence_filter_order) = state.namelists.turbulence
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid

    if nbz < turbulence_filter_order
        error("Error in smooth_tke!: nbz < turbulence_filter_order!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        output[i, j, k] =
            sum(
                input[i, j, kk] * jac[i, j, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order)
            ) / sum(
                jac[i, j, kk] for kk in
                (k - turbulence_filter_order):(k + turbulence_filter_order)
            )
    end

    return
end
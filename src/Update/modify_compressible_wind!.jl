function modify_compressible_wind!(state::State, operation::Function)
    (; model) = state.namelists.setting
    modify_compressible_wind!(state, operation, model)
    return
end

function modify_compressible_wind!(
    state::State,
    operation::Function,
    model::AbstractModel,
)
    return
end

function modify_compressible_wind!(
    state::State,
    operation::Function,
    model::Compressible,
)
    (; namelists, domain) = state
    (; nbz) = namelists.domain
    (; zboundaries) = namelists.setting
    (; nxx, nyy, i0, i1, j0, j1, k0, k1) = domain
    (; jac) = state.grid
    (; u, v, w, p) = state.variables.predictands

    if zboundaries != SolidWallBoundaries()
        error("Error in modify_compressible_wind!: Unknown zboundaries!")
    end

    # Modify the zonal wind.
    for kz in k0:k1, jy in 1:nyy, ix in i0:i1
        u[ix, jy, kz] = operation(
            u[ix, jy, kz],
            (
                jac[ix, jy, kz] * p[ix, jy, kz] +
                jac[ix + 1, jy, kz] * p[ix + 1, jy, kz]
            ) / 2,
        )
    end
    set_zonal_boundaries_of_field!(u, namelists, domain)
    for k in 1:nbz
        @views u[:, :, k0 - k] .= u[:, :, k0 + k - 1]
        @views u[:, :, k1 + k] .= u[:, :, k1 - k + 1]
    end

    # Modify the meridional wind.
    for kz in k0:k1, jy in j0:j1, ix in 1:nxx
        v[ix, jy, kz] = operation(
            v[ix, jy, kz],
            (
                jac[ix, jy, kz] * p[ix, jy, kz] +
                jac[ix, jy + 1, kz] * p[ix, jy + 1, kz]
            ) / 2,
        )
    end
    set_meridional_boundaries_of_field!(v, namelists, domain)
    for k in 1:nbz
        @views v[:, :, k0 - k] .= v[:, :, k0 + k - 1]
        @views v[:, :, k1 + k] .= v[:, :, k1 - k + 1]
    end

    # Modify the transformed vertical wind.
    for kz in k0:k1, jy in 1:nyy, ix in 1:nxx
        w[ix, jy, kz] = operation(
            w[ix, jy, kz],
            jac[ix, jy, kz] *
            jac[ix, jy, kz + 1] *
            (p[ix, jy, kz] + p[ix, jy, kz + 1]) /
            (jac[ix, jy, kz] + jac[ix, jy, kz + 1]),
        )
    end
    w[:, :, k0 - 1] .= 0.0
    w[:, :, k1] .= 0.0
    for k in 1:nbz
        @views w[:, :, k0 - k] .= -w[:, :, k0 + k - 2]
        @views w[:, :, k1 + k] .= -w[:, :, k1 - k]
    end

    return
end
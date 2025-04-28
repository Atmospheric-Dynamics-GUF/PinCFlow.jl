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
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; u, v, w, p) = state.variables.predictands

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        u[ix, jy, kz] = operation(
            u[ix, jy, kz],
            (
                jac[ix, jy, kz] * p[ix, jy, kz] +
                jac[ix + 1, jy, kz] * p[ix + 1, jy, kz]
            ) / 2,
        )

        v[ix, jy, kz] = operation(
            v[ix, jy, kz],
            (
                jac[ix, jy, kz] * p[ix, jy, kz] +
                jac[ix, jy + 1, kz] * p[ix, jy + 1, kz]
            ) / 2,
        )

        w[ix, jy, kz] = operation(
            w[ix, jy, kz],
            jac[ix, jy, kz] *
            jac[ix, jy, kz + 1] *
            (p[ix, jy, kz] + p[ix, jy, kz + 1]) /
            (jac[ix, jy, kz] + jac[ix, jy, kz + 1]),
        )
    end

    return
end
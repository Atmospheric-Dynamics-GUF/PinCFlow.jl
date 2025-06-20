"""
    modify_compressible_wind!(state::State, operation::Function)

Apply operation to wind fields based on model type.

This function dispatches to the appropriate model-specific implementation
for modifying wind components in compressible flow models.

# Arguments

  - `state::State`: Simulation state containing model configuration
  - `operation::Function`: Binary operation to apply (typically `*` or `/`)
"""
function modify_compressible_wind!(state::State, operation::Function)
    (; model) = state.namelists.setting
    modify_compressible_wind!(state, operation, model)
    return
end

"""
    modify_compressible_wind!(state::State, operation::Function, model::AbstractModel)

No-op for non-compressible models.

Non-compressible models don't require wind field modifications as they
don't use pressure-density coupling in the momentum equations.

# Arguments

  - `state::State`: Simulation state (unused)
  - `operation::Function`: Operation function (unused)
  - `model::AbstractModel`: Non-compressible model type
"""
function modify_compressible_wind!(
    state::State,
    operation::Function,
    model::AbstractModel,
)
    return
end

"""
    modify_compressible_wind!(state::State, operation::Function, model::Compressible)

Modify wind components for compressible model using pressure weighting.

Applies the specified operation to wind components using pressure-weighted
averages. This is used in semi-implicit time stepping to account for
pressure-density coupling in compressible flows.

# Arguments

  - `state::State`: Simulation state containing grid, domain, and predictand fields
  - `operation::Function`: Binary operation (typically `*` for scaling, `/` for unscaling)
  - `model::Compressible`: Compressible model type

# Details

  - Modifies u-component using pressure average between adjacent x-points
  - Modifies v-component using pressure average between adjacent y-points
  - Modifies w-component using Jacobian-weighted pressure average between z-levels
  - Operations are applied only within the computational domain (i0:i1, j0:j1, k0:k1)
"""
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

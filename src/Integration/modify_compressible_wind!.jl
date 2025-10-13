"""
```julia
modify_compressible_wind!(state::State, operation::Function)
```

Modify the wind with ``J P`` if the atmosphere is compressible by dispatching to the appropriate method.

```julia
modify_compressible_wind!(
    state::State,
    operation::Function,
    model::Union{Boussinesq, PseudoIncompressible},
)
```

Return in non-compressible modes.

```julia
modify_compressible_wind!(
    state::State,
    operation::Function,
    model::Compressible,
)
```

Interpolate ``J P`` to the wind grids and replace the wind components with the result of applying `operation` to them and the interpolations.

# Arguments

  - `state`: Model state.

  - `operation`: Binary operation used for modification.

  - `model`: Dynamic equations.
"""
function modify_compressible_wind! end

function modify_compressible_wind!(state::State, operation::Function)
    (; model) = state.namelists.setting
    modify_compressible_wind!(state, operation, model)
    return
end

function modify_compressible_wind!(
    state::State,
    operation::Function,
    model::Union{Boussinesq, PseudoIncompressible},
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

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        u[i, j, k] = operation(
            u[i, j, k],
            (jac[i, j, k] * p[i, j, k] + jac[i + 1, j, k] * p[i + 1, j, k]) / 2,
        )

        v[i, j, k] = operation(
            v[i, j, k],
            (jac[i, j, k] * p[i, j, k] + jac[i, j + 1, k] * p[i, j + 1, k]) / 2,
        )

        w[i, j, k] = operation(
            w[i, j, k],
            jac[i, j, k] * jac[i, j, k + 1] * (p[i, j, k] + p[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1]),
        )
    end

    return
end

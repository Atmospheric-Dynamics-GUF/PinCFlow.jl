function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_wind_factor(state, indices, variable, model)
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    model::AbstractModel,
)
    return 1.0
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix + 1, jy, kz] * p[ix + 1, jy, kz]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix, jy + 1, kz] * p[ix, jy + 1, kz]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return jac[ix, jy, kz] *
           jac[ix, jy, kz + 1] *
           (p[ix, jy, kz] + p[ix, jy, kz + 1]) /
           (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
end

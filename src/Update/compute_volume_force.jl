function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; testcase) = state.namelists.setting
    return compute_volume_force(state, indices, variable, testcase)
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)
    return 0.0
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    testcase::AbstractWKBTestCase,
)
    (; dudt) = state.wkb.tendencies
    (ix, jy, kz) = indices
    return (dudt[ix, jy, kz] + dudt[ix + 1, jy, kz]) / 2
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)
    (; dvdt) = state.wkb.tendencies
    (ix, jy, kz) = indices
    return (dvdt[ix, jy, kz] + dvdt[ix, jy + 1, kz]) / 2
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    testcase::AbstractWKBTestCase,
)
    (; met) = state.grid
    (; dudt, dvdt) = state.wkb.tendencies
    (; ix, jy, kz) = indices
    return (
        met[ix, jy, kz, 1, 3] * dudt[ix, jy, kz] +
        met[ix, jy, kz, 2, 3] * dvdt[ix, jy, kz] +
        met[ix, jy, kz + 1, 1, 3] * dudt[ix, jy, kz + 1] +
        met[ix, jy, kz + 1, 2, 3] * dvdt[ix, jy, kz + 1]
    ) / 2
end
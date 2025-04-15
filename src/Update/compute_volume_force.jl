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
    return state.wkb.integrals.dudt[indices...]
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)
    return state.wkb.integrals.dvdt[indices...]
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    testcase::AbstractWKBTestCase,
)
    return state.wkb.integrals.dwdt[indices...]
end
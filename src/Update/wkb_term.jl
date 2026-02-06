function wkb_term end

function wkb_term(state::State, i::Integer, j::Integer, k::Integer)
    (; wkb_mode) = state.namelists.wkb

    return wkb_term(state, i, j, k, wkb_mode)
end

function wkb_term(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    wkb_mode::NoWKB,
)
    return 0.0
end

function wkb_term(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)
    (; dtkedt) = state.wkb.tendencies

    return dtkedt[i, j, k]
end

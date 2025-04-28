function set_compressible_meridional_boundaries!(
    state::State,
    model::AbstractModel,
)
    return
end

function set_compressible_meridional_boundaries!(
    state::State,
    model::Compressible,
)
    (; namelists, domain) = state
    (; p) = state.variables.predictands
    set_meridional_boundaries_of_field!(p, namelists, domain)
    return
end
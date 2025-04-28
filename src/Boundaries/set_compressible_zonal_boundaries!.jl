function set_compressible_zonal_boundaries!(
    state::State,
    model::AbstractModel,
)
    return
end

function set_compressible_zonal_boundaries!(
    state::State,
    model::Compressible,
)
    (; namelists, domain) = state
    (; p) = state.variables.predictands
    set_zonal_boundaries_of_field!(p, namelists, domain)
    return
end
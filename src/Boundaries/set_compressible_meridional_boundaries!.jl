"""
    set_compressible_meridional_boundaries!(state, model::AbstractModel)

No-op for non-compressible models.
"""
function set_compressible_meridional_boundaries!(
    state::State,
    model::AbstractModel,
)
    return
end

"""
    set_compressible_meridional_boundaries!(state, model::Compressible)

Set meridional boundaries for pressure field in compressible model.
"""
function set_compressible_meridional_boundaries!(
    state::State,
    model::Compressible,
)
    (; namelists, domain) = state
    (; p) = state.variables.predictands
    set_meridional_boundaries_of_field!(p, namelists, domain)
    return
end
function save_backups!(state::State, variables::Vararg{Symbol})
    (; backups, predictands) = state.variables

    for field in variables
        getfield(backups, Symbol(field, :old)) .= getfield(predictands, field)
    end

    return
end

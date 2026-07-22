# A singleton struct that holds the names of the active fields in its type parameters depending on the chosen nucleation type.

struct IceActiveVars{Fields}
    # Empty singleton struct
end

"""
The function `get_IceActiveVars(icepredictands)` checks which arrays inside `icepredictands` 
are not empty and returns an `IceActiveVars` struct parameterized with a tuple of those active field names.
"""

function get_IceActiveVars(icepredictands)
    # Find all field names (symbols) where the corresponding array is not empty.
    active_symbols = [f for f in fieldnames(typeof(icepredictands)) if !isempty(getfield(icepredictands, f))]

    # Return the struct with the active symbols as its type parameter 
    return IceActiveVars{Tuple(active_symbols)}()
end

function ice_active_vars_tuple(::IceActiveVars{Fields}) where {Fields}
    return Fields
end
"""
```julia
dispatch(values::Expr, input::Any)::Any
```

If `input` is an expression, find the first dynamic value dispatch in it and make all similar dynamic value dispatches static by inserting `input` into an `if` block which compares the found argument to the options in `values` and replacing it accordingly in all `Val` calls.

# Arguments

  - `values`: Expression of a tuple of allowed values.

  - `input`: Input object.

# See also

  - [`PinCFlow.Macros.find_argument`](@ref)

  - [`PinCFlow.Macros.replace_argument`](@ref)
"""
function dispatch end

function dispatch(values::Expr, input::Any)::Any
    if values.head !== :tuple
        error("Tuple expected for values!")
    end

    argument = find_argument(input)

    if argument !== nothing
        condition = ""
        for (index, value) in enumerate(Core.eval(@__MODULE__, values))
            if value isa Symbol
                literal = ":$(value)"
            elseif value isa AbstractString
                literal = "\"$(value)\""
            else
                literal = "$(value)"
            end
            prefix = index == 1 ? "if" : "elseif"
            condition *=
                "$(prefix) $(argument) === $(literal)\n" *
                string(replace_argument(input, argument, Meta.parse(literal))) *
                "\n"
        end
        condition *= "else\nerror(\"Invalid $(argument) option!\")\nend"

        return Meta.parse(condition)
    else
        return input
    end
end

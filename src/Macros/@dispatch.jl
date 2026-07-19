"""
```julia
@dispatch(values::Expr, input::Any)
```

Find the first dynamic value dispatch in `input` and make all dynamic value dispatches with the same argument static.

The scope in which this macro is to be applied must include the argument of the first dynamic value dispatch, and the available methods of the involved functions must be covered by the options given in `values`.

# Arguments

  - `values`: Expression of a tuple of allowed values.

  - `input`: Input object.

# See also

  - [`PinCFlow.Macros.dispatch`](@ref)
"""
macro dispatch end

macro dispatch(values::Expr, input::Any)
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

        return esc(Meta.parse(condition))
    else
        return esc(input)
    end
end

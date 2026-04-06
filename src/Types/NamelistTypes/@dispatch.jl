"""
```julia
@dispatch(values::Expr, input::Expr)
```

Macro that makes value dispatch static for a the given list of values.

# Arguments

  - `values`: Expression of a tuple of allowed values.

  - `input`: Input expression with `Val` calls.
"""
macro dispatch end

macro dispatch(values::Expr, input::Expr)
    code = string(input)
    range = findfirst(r"\bVal\((?!\s*(:\w+\b|\d+\b|true\b|false\b))", code)

    @ivy if range !== nothing
        start = first(range) + 4
        stop = Meta.parse(code, start - 1; greedy = false)[2] - 2
        parameter = code[start:stop]

        startswith(parameter, r"\s*\W") &&
            error("Invalid Val call for @dispatch!")

        modified_code = ""
        for (index, value) in enumerate(Core.eval(@__MODULE__, values))
            if typeof(value) <: Symbol
                literal = ":$(value)"
            elseif typeof(value) <: AbstractString
                literal = "\"$(value)\""
            else
                literal = "$(value)"
            end
            prefix = index == 1 ? "if" : "elseif"
            modified_code *=
                "$(prefix) $(parameter) === $(literal)\n" *
                replace(
                    code,
                    Regex(
                        "\\bVal\\(\\s*$(parameter)\\s*\\)",
                    ) => "Val($(literal))",
                ) *
                "\n"
        end
        modified_code *= "else\nerror(\"Invalid $(parameter) option!\")\nend"
    else
        modified_code = code
    end

    return esc(Meta.parse(modified_code))
end

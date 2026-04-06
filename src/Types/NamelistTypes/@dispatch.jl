"""
```julia
@dispatch(input::Expr)
```

Macro that makes value dispatch static for a predefined list of symbols.

# Symbols

  - `:monotone_centered_variant`

# Arguments

  - `input`: Input expression with a `Val` call.
"""
macro dispatch end

macro dispatch(input::Expr)
    symbols = (:monotone_centered_variant,)

    code = string(input)
    calls = match(r"\bVal\(", code)

    modified_code = ""
    if calls !== nothing
        length(calls.offsets) > 1 &&
            error("@dispatch is not ready for multiple Val calls yet!")

        start = calls.offset + 4
        stop = Meta.parse(code, start - 1; greedy = false)[2] - 2
        parameter = code[start:stop]

        for (index, symbol) in enumerate(symbols)
            index > 1 && (modified_code *= "else")
            modified_code *=
                "if $(parameter) == :$(symbol)\n" *
                code[1:(start - 1)] *
                ":$(symbol)" *
                code[(stop + 1):end] *
                "\n"
        end
        modified_code *= "else\nerror(\"Unknown option!\")\nend"
    end

    return esc(Meta.parse(modified_code))
end

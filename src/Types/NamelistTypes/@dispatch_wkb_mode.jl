"""
```julia
@dispatch_wkb_mode(input::Expr)
```

Macro that makes value dispatch static for the `wkb_mode` parameter of `WKBNamelist`.

The parameter can take any of the following values:

  - `:NoWKB`

  - `:SteadyState`

  - `:SingleColumn`

  - `:MultiColumn`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_wkb_mode end

macro dispatch_wkb_mode(input::Expr)
    return esc(
        quote
            @dispatch (:NoWKB, :SteadyState, :SingleColumn, :MultiColumn) $(
                input
            )
        end,
    )
end

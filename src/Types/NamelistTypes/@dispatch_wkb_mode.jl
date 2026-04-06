"""
```julia
@dispatch_wkb_mode(input::Expr)
```

Macro that makes value dispatch for the `wkb_mode` parameter static.

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

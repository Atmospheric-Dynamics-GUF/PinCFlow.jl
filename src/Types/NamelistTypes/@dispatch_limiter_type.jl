"""
```julia
@dispatch_limiter_type(input::Expr)
```

Macro that makes value dispatch static for the `limiter_type` parameter of `DiscretizationNamelist`.

The parameter can take any of the following values:

  - `:MCVariant`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_limiter_type end

macro dispatch_limiter_type(input::Expr)
    return esc(quote
        @dispatch (:MCVariant,) $(input)
    end)
end

"""
```julia
@dispatch_limiter_type(input::Expr)
```

Macro that makes value dispatch for the `limiter_type` parameter static.

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_limiter_type end

macro dispatch_limiter_type(input::Expr)
    return esc(quote
        @dispatch (:monotone_centered_variant,) $(input)
    end)
end

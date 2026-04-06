"""
```julia
@dispatch_filter_order(input::Expr)
```

Macro that makes value dispatch for the `filter_order` parameter static.

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_filter_order end

macro dispatch_filter_order(input::Expr)
    return esc(quote
        @dispatch (1, 2, 3, 4) $(input)
    end)
end

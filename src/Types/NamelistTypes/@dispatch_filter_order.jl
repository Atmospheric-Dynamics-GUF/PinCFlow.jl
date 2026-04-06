"""
```julia
@dispatch_filter_order(input::Expr)
```

Macro that makes value dispatch static for the `filter_order` parameter of `WKBNamelist`.

The parameter can take any of the following values:

  - `1`

  - `2`

  - `3`

  - `4`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_filter_order end

macro dispatch_filter_order(input::Expr)
    return esc(quote
        @dispatch (1, 2, 3, 4) $(input)
    end)
end

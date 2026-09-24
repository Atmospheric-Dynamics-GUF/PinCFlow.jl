"""
```julia
@dispatch_filter_type(input::Expr)
```

Macro that makes value dispatch static for the `filter_type` parameter of `WKBNamelist`.

The parameter can take any of the following values:

  - `:BoxFilter`

  - `:ShapiroFilter`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_filter_type end

macro dispatch_filter_type(input::Expr)
    return esc(quote
        @dispatch (:BoxFilter, :ShapiroFilter) $(input)
    end)
end

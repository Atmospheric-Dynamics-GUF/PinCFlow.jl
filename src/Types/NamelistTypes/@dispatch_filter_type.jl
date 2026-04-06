"""
```julia
@dispatch_filter_type(input::Expr)
```

Macro that makes value dispatch for the `filter_type` parameter static.

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_filter_type end

macro dispatch_filter_type(input::Expr)
    return esc(quote
        @dispatch (:BoxFilter, :ShapiroFilter) $(input)
    end)
end

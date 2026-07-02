"""
```julia
@dispatch_orographic_flow(input::Expr)
```

Macro that makes value dispatch static for the `orographic_flow` parameter of `WKBNamelist`.

The parameter can take any of the following values:

  - `:Surface`

  - `:Summit`

  - `:Average`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_background end

macro dispatch_background(input::Expr)
    return esc(quote
        @dispatch (:Surface, :Summit, :Average) $(input)
    end)
end

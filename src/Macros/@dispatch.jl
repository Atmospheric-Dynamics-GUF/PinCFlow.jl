"""
```julia
@dispatch(values::Expr, input::Expr)
```

Find the first dynamic value dispatch in `input` and make all similar dynamic value dispatches static.

The scope in which this macro is to be applied must include the argument of the first dynamic value dispatch.

# Arguments

  - `values`: Expression of a tuple of allowed values.

  - `input`: Input expression with dynamic value dispatches.
"""
macro dispatch end

macro dispatch(values::Expr, input::Expr)
    return esc(dispatch(values, input))
end

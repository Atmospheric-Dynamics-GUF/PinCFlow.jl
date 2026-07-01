"""
```julia
@dispatch(values::Expr, input::Expr)
```

Find the first dynamic value dispatch in `input` and make all dynamic value dispatches with the same argument static.

The scope in which this macro is to be applied must include the argument of the first dynamic value dispatch, and the available methods of the involved functions must be covered by the options given in `values`.

# Arguments

  - `values`: Expression of a tuple of allowed values.

  - `input`: Input expression with dynamic value dispatches.

# See also

  - [`PinCFlow.Macros.dispatch`](@ref)
"""
macro dispatch end

macro dispatch(values::Expr, input::Expr)
    return esc(dispatch(values, input))
end

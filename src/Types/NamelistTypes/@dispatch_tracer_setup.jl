"""
```julia
@dispatch_tracer_setup(input::Expr)
```

Macro that makes value dispatch static for the `tracer_setup` parameter of `TracerNamelist`.

The parameter can take any of the following values:

  - `:NoTracer`

  - `:TracerOn`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_tracer_setup end

macro dispatch_tracer_setup(input::Expr)
    return esc(quote
        @dispatch (:NoTracer, :TracerOn) $(input)
    end)
end

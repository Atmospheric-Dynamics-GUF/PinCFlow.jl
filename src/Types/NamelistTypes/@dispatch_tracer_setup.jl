"""
```julia
@dispatch_tracer_setup(input::Expr)
```

Macro that makes value dispatch for the `tracer_setup` parameter static.

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_tracer_setup end

macro dispatch_tracer_setup(input::Expr)
    return esc(quote
        @dispatch (:no_tracer, :tracer_on) $(input)
    end)
end

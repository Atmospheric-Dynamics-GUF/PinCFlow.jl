"""
```julia
@dispatch_model(input::Expr)
```

Macro that makes value dispatch for the `model` parameter static.

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_model end

macro dispatch_model(input::Expr)
    return esc(
        quote
            @dispatch (:Boussinesq, :PseudoIncompressible, :Compressible) $(
                input
            )
        end,
    )
end

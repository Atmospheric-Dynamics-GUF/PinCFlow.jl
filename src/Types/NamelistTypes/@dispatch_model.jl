"""
```julia
@dispatch_model(input::Expr)
```

Macro that makes value dispatch static for the `model` parameter of `AtmosphereNamelist`.

The parameter can take any of the following values:

  - `:Boussinesq`

  - `:PseudoIncompressible`

  - `:Compressible`

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

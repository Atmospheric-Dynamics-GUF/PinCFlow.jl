"""
```julia
@dispatch_turbulence_scheme(input::Expr)
```

Macro that makes value dispatch static for the `turbulence_scheme` parameter of `TurbulenceNamelist`.

The parameter can take any of the following values:

  - `:NoTurbulence`

  - `:TKEScheme`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_turbulence_scheme end

macro dispatch_turbulence_scheme(input::Expr)
    return esc(quote
        @dispatch (:NoTurbulence, :TKEScheme) $(input)
    end)
end

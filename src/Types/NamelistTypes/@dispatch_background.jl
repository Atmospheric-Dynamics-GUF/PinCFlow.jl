"""
```julia
@dispatch_background(input::Expr)
```

Macro that makes value dispatch static for the `background` parameter of `AtmosphereNamelist`.

The parameter can take any of the following values:

  - `:NeutralStratification`

  - `:StableStratification`

  - `:Isothermal`

  - `:Isentropic`

  - `:Realistic`

  - `:LapseRates`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_background end

macro dispatch_background(input::Expr)
    return esc(
        quote
            @dispatch (
                :NeutralStratification,
                :StableStratification,
                :Isothermal,
                :Isentropic,
                :Realistic,
                :LapseRates,
            ) $(input)
        end,
    )
end

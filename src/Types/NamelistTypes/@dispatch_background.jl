"""
```julia
@dispatch_background(input::Expr)
```

Macro that makes value dispatch for the `background` parameter static.

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

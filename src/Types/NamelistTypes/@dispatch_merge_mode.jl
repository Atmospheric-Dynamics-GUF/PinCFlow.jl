"""
```julia
@dispatch_merge_mode(input::Expr)
```

Macro that makes value dispatch static for the `merge_mode` parameter of `WKBNamelist`.

The parameter can take any of the following values:

  - `:ConstantWaveAction`

  - `:ConstantWaveEnergy`

# Arguments

  - `input`: Input expression with `Val` calls.
"""
macro dispatch_merge_mode end

macro dispatch_merge_mode(input::Expr)
    return esc(
        quote
            @dispatch (:ConstantWaveAction, :ConstantWaveEnergy) $(input)
        end,
    )
end

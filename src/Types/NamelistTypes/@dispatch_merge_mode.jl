"""
```julia
@dispatch_merge_mode(input::Expr)
```

Macro that makes value dispatch for the `merge_mode` parameter static.

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

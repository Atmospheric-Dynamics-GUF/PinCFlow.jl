"""
```julia
reset_thomas!(state::State)
```

Reset arrays needed for the Thomas tridiagonal algorithm to ``0``.

# Arguments

  - `state`: Model state.
"""
function reset_thomas! end

function reset_thomas!(state::State)
    (; auxiliaries) = state.variables

    for field_name in (:ath, :bth, :cth, :fth, :qth, :qth_bc, :fth_bc)
        field = getfield(auxiliaries, field_name)
        @share field = 0.0
    end

    return
end

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

    for field in (:ath, :bth, :cth, :fth, :qth, :pth, :qth_bc, :fth_bc)
        getfield(auxiliaries, field) .= 0.0
    end

    return
end

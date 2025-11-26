function reset_thomas! end

function reset_thomas!(state::State)
    (; auxiliaries) = state.variables

    for field in (
        :athglob,
        :bthglob,
        :cthglob,
        :fthglob,
        :qthglob,
        :pthglob,
        :qthglob_bc,
        :fthglob_bc,
    )
        getfield(auxiliaries, field) .= 0.0
    end

    return
end

function reset_thomas! end

function reset_thomas!(state::State)

    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

    ath .= 0.0
    bth .= 0.0
    cth .= 0.0
    fth .= 0.0
    qth .= 0.0
    pth .= 0.0
    qth_bc .= 0.0
    fth_bc .= 0.0

    return 
end
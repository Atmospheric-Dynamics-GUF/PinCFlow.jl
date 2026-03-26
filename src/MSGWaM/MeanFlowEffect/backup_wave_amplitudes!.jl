function backup_wave_amplitudes! end

function backup_wave_amplitudes!(state::State)
    (; uhatold, vhatold, whatold, bhatold, pihatold, chihatold) = state.wkb.wkbauxiliaries
    (; uhat, vhat, what, bhat, pihat, chihat) = state.wkb.integrals

    uhatold .= uhat
    vhatold .= vhat
    whatold .= what
    bhatold .= bhat
    pihatold .= pihat
    chihatold .= chihat

    return
end

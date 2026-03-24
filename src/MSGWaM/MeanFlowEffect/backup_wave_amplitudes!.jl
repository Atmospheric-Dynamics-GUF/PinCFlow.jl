function backup_wave_amplitudes! end

function backup_wave_amplitudes!(state::State)
    (; uhatold, vhatold, whatold, bhatold, pihatold) = state.wkb.wkbauxiliaries
    (; uhat, vhat, what, bhat, pihat) = state.wkb.integrals

    uhatold .= uhat
    vhatold .= vhat
    whatold .= what
    bhatold .= bhat
    pihatold .= pihat

    return
end
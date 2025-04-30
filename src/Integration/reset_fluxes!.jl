function reset_fluxes!(state::State, fluxes::Fluxes)
    (; phirho, phirhop, phiu, phiv, phiw) = state.variables.fluxes

    phirho .= fluxes.phirho
    phirhop .= fluxes.phirhop
    phiu .= fluxes.phiu
    phiv .= fluxes.phiv
    phiw .= fluxes.phiw

    return
end

function check_tke!(state::State)

    (; k0, k1, j0, j1, i0, i1) = state.domain
    (; tke) = state.turbulence.turbulencepredictands 
    (; tkemin) = state.turbulence.turbulenceconstants
    (; rhobar) = state.atmosphere

    for k in k0:k1, j in j0:j1, i in i0:i1 
        tke[i, j, k] = max(tke[i, j, k], tkemin * rhobar[i, j, k])
    end

    return

end
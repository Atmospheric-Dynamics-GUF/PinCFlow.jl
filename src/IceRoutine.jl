function compute_source_ice!(state::State)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    
    (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands
    (; iceconstants) = state.ice
    (; icesource) = state.ice
    (; iceauxiliaries) = state.ice
    (; kappainv, pref, gamma) = state.constants
    (; press0_dim) = state.namelists.atmosphere

    p0 = press0_dim / pref

    for k in k0:k1, j in j0:j1, i in i0:i1
        
        # Question exn_p = pi(i, j, k) + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
        exn_p = pip[i, j, k] + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)

        rqv = state.ice.icepredictands.qv[i, j, k] 
        pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
        rhoMean = rhostrattfc[i, j, k] 
        rho_full = rho[i, j, k] + rhoMean
        theta = pstrattfc[i, j, k] / rho_full

        temp = theta * exn_p

        psi = psat_ice(temp, iceconstants)

        NIce = state.ice.icepredictands.n[i, j, k] # N_v = \rho n

        sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

        iceauxiliaries.iaux1[i, j, k] = sice
        iceauxiliaries.iaux2[i, j, k] = psi
        iceauxiliaries.iaux3[i, j, k] = rqv

        if sice >= iceconstants.S_c  
            icesource.nsource[i, j, k] = dot_n(sice, rho_full, iceconstants)
        else
            icesource.nsource[i, j, k] = 0.0       
        end 

        dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)

        icesource.qvsource[i, j, k] = dqv
        icesource.qsource[i, j, k] = -dqv
    end 

    return
end


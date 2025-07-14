module IceRoutine

using ..Types
#using ..IceConstants

function compute_source_ice!(state::State)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    
    (; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
    (; rho, rhop, u, v, w, pip, p) = state.variables.predictands
    (; iceconstants) = state.ice
    (; icesource) = state.ice
    (; kappainv) = state.constants

    for k in k0:k1, j in j0:j1, i in i0:i1
    
        rqv = state.ice.predictands.qv[i, j, k] 
        pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
        rhoMean = rhostrattfc[i, j, k] 
        rho_full = rho[i, j, k] + rhoMean
        
        sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

        if sice >= iceconstants.S_c  
            icesource.nsource[i, j, k] = dot_n(sice, rho_full, iceconstants)
        else
            icesource.nsource[i, j, k] = 0.0       
        end 

        dqv = dot_qv(SIce, NIce, temp, pres, psi, iceconstants)

        icesource.qvsource[i, j, k] = dqv
        icesource.qvsource[i, j, k] = -dqv
    end 

    return
end

function sat_ratio(Qv::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, rhoMean::AbstractFloat, cons::IceConstants)
    (; epsil0hat) = cons
    return Qv / rhoMean * p / p_si / epsil0hat
end     

function dot_qv(S::AbstractFloat, N::AbstractFloat, T::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, cons::IceConstants)
    (; Dep) = cons
    return - Dep * (S - 1.0) * T * N * p_si / p
end    

function dot_n(S::AbstractFloat, rho::AbstractFloat, cons::IceConstants)
    (; J_nuc, B_nuc, S_c) = cons
    return rho * J_nuc * exp(B_nuc * (S - S_c))
end

export compute_source_ice!, sat_ratio, dot_qv, dot_n

end
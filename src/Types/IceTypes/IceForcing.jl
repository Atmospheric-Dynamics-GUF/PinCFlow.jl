mutable struct IceForcing{A <: AbstractFloat}
    #tau_q_sink :: A # time scale for q sink (s)
    #tau_n_sink :: A # time scale for n sink (s)
    #tau_qv_source :: A # time scale for qv source (s)

    time_physical :: A # physical time (s)
end

function IceForcing(constants::Constants)
    (; tref, rhoref, lref) =  constants

    tRef = tref
    rhoRef = rhoref
    lRef = lref

    mRef = rhoRef * lRef ^ 3     # reference mass

    #tau_q_sink = 3.0e-11 * tRef * mRef^(2.0/3.0)
    #tau_n_sink = 2.0 * tau_q_sink
    # tau_qv_source = 3.0e-11 * tRef to be added later

    time_physical = 0.0 # to be updated during the simulation

    return IceForcing(time_physical)
end
mutable struct IceForcing{
    A <: AbstractFloat,
    B <: AbstractVector{<:AbstractFloat}
}
    time_physical :: A    # actual physical time in seconds
    qv_ref :: B           # Vertical reference profile of qv
end


function IceForcing(
    constants::Constants,
    domain::Domain
)

    (; tref, rhoref, lref) =  constants
    (; nzz) = domain

    tRef = tref
    rhoRef = rhoref
    lRef = lref

    mRef = rhoRef * lRef ^ 3     # reference mass

    #tau_q_sink = 3.0e-11 * tRef * mRef^(2.0/3.0)
    #tau_n_sink = 2.0 * tau_q_sink
    # tau_qv_source = 3.0e-11 * tRef to be added later

    time_physical = 0.0 # to be updated during the simulation
    qv_ref = zeros(Float64, nzz)

    return IceForcing(time_physical, qv_ref)
end
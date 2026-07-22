mutable struct IceForcing{
    A <: AbstractFloat,
    B <: AbstractVector{<:AbstractFloat}
}
    time_physical :: A    # actual physical time in seconds
    qv_ref :: B           # Vertical reference profile of qv
end


function IceForcing(
    constants::Constants, # constants can be removed
    domain::Domain
)

    (; nzz) = domain

    time_physical = 0.0 # physical time in seconds
    qv_ref = zeros(Float64, nzz)

    return IceForcing(time_physical, qv_ref)
end
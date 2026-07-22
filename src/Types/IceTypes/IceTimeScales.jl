"""
    IceTimeScales.jl

    This file contains the sedimentation and relaxation timescales used in the ice microphysics model.
"""
struct IceTimeScales{A <: AbstractFloat}
    
    tau_sink::A
    tau_relax::A

end    

function IceTimeScales()

    tau_sink = 3.3e2 # s
    tau_relax = 3.0e3 # s

    return IceTimeScales{Float64}(
        tau_sink,
        tau_relax
    )
end
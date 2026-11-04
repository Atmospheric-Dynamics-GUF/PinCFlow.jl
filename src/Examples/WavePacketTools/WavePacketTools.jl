"""
```julia
WavePacketTools
```

Module that defines helper functions for setting up the wave-packet and WKB-wave-packet example simulations.

# Dependencies

  - [`PinCFlow.Macros`](@ref)

  - [`PinCFlow.Types`](@ref)

# Exports

  - [`PinCFlow.Examples.WavePacketTools.bhat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.n2`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.omega`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.phi`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.pihat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.rhobar`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.thetabar`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.uhat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.vhat`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.wave_action_density`](@ref)

  - [`PinCFlow.Examples.WavePacketTools.what`](@ref)
"""
module WavePacketTools

using ...Macros
using ...Types

include("bhat.jl")
include("envelope.jl")
include("ijk.jl")
include("n2.jl")
include("omega.jl")
include("phi.jl")
include("pihat.jl")
include("rhobar.jl")
include("thetabar.jl")
include("uhat.jl")
include("vhat.jl")
include("wave_action_density.jl")
include("what.jl")

export bhat,
    n2,
    omega,
    phi,
    pihat,
    rhobar,
    thetabar,
    uhat,
    vhat,
    wave_action_density,
    what

end

"""
```julia
TriadInteractions
```

Module for computing the scattering integral of gravity waves.

Provides functions that compute scattering integral by integrating the RHS of kinetic equation in the 2D spectral space (kp, m) by assuming the horizontal isotropy and vertical symmetry N(kp, |m|). 

  - [`PinCFlow.Types`](@ref)

  - [`PinCFlow.Boundaries`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation`](@ref)

  - [`PinCFlow.MSGWaM.RayUpdate`](@ref)
"""
module TriadInteractions

using ..MeanFlowEffect
using ..Interpolation
using ..RayUpdate
using ...Types
using ...Boundaries
using ...PinCFlow



include("compute_spectral_cell_indices.jl")
include("compute_wave_spectrum!.jl")
include("half_logwidth.jl")



export compute_wave_spectrum!

end

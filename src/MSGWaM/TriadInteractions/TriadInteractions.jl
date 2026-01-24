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
include("half_logwidth.jl")
include("get_wave_spectrum!.jl")
include("apply_triad_interactions!.jl")
include("interaction_matrix.jl")
include("compute_delta_pq.jl")
include("compute_delta.jl")
include("compute_g_prime.jl")
include("compute_kp1kp2.jl")
include("compute_m1m2.jl")
include("trapazoidal_with_logbin.jl")
include("update_interpolation_coef!.jl")
include("interpolate_nk.jl")
include("update_wave_spectrum!.jl")
include("compute_st_k.jl")
include("compute_scattering_integral!.jl")
include("get_ray_volumes!.jl")
include("launch_new_ray_vol!.jl")
include("initialize_wave_spectrum!.jl")
include("check_resonance.jl")

export get_wave_spectrum!,
       apply_triad_interactions!,
       initialize_wave_spectrum!

end

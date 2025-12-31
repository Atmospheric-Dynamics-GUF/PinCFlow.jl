"""
```julia
TriadTendencies{A <: AbstractVector{<:AbstractFloat}, B <: AbstractArray{<: AbstractFloat, 5}}
```

Triad interaction fields.

```julia
TriadTendencies{A <: AbstractVector{<:AbstractFloat}, B <: AbstractArray{<: AbstractFloat, 5}}
    kp::A
    m::A
    wavespectrum::B
    st_k::B
    col_int::B
    kin_box::KinematicBox
    interp_coef::InterpCoef
```

Construct a `TriadTendencies` instance, with arrays sized according to the given dimensions.

# Fields

  - `kp::A`: horizontal wave number grid.

  - `m::A`: vertical wave number grid.

  - `wavespectrum::B`: wave action spectrum in 5D phase space.
  
  - `st_k::B`: integrand of collision integral.

  - `col_int::B`: collision integral.

  - `kin_box::KinematicBox`: kinematic box container.


# Arguments

  - `domain`

  - `constants`
  
  - `wkb_mode`
  
"""

const RaySignature = Tuple{Int, Int, Int, Int, Float64, Float64}


struct TriadTendencies{A <: AbstractArray{<: AbstractFloat, 5}, B <: AbstractArray{<: AbstractFloat, 2}}
    spec_grid::SpectralGrid
    wavespectrum::A
    ray_vol_signature::Array{Vector{RaySignature},5}
    col_int::B
    kin_box::KinematicBox
    interp_coef::InterpCoef
end


function TriadTendencies(namelists::Namelists,
   domain::Domain,
   wkb_mode::NoWKB)::TriadTendencies

  spec_grid = SpectralGrid(wkb_mode)
   
  kin_box = KinematicBox(wkb_mode)

  interp_coef = InterpCoef(wkb_mode)  
  
  v = Array{Vector{RaySignature}}(undef, 0, 0, 0, 0, 0)
  
  return TriadTendencies(spec_grid, zeros(0, 0, 0, 0, 0), v, zeros(0, 0),
    kin_box, interp_coef)
  end

function TriadTendencies(namelists::Namelists,
   domain::Domain,
   constants::Constants,
   nray_max::Integer,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})::TriadTendencies



   (; nxx, nyy, nzz) = domain

   # Compute the grid in kp and m direction
  

   spec_grid = SpectralGrid(namelists, constants, wkb_mode)
   
   (; kp, m) = spec_grid

   println(kp, m)
   kpl = length(kp)
   ml = length(m)
   kpmin = kp[1] 
   kpmax = kp[end]
  
   amin = (kpmin / kpl) .* ones(kpl)
   amax = Float64.(kp)
   ma = Int.(max.(8*ones(kpl), 1:kpl))  
  
   mq = 2 * kpl
   qmin = kpmin / mq
   qmax = 2.0 * kpmax
   kin_box = KinematicBox(amin, amax, ma, qmin, qmax, mq, wkb_mode)

   interp_coef = InterpCoef(kp, m, wkb_mode)

   ray_vol_signature = Array{Vector{RaySignature},5}(undef, nxx, nyy, nzz, kpl, ml)

  for idx in eachindex(ray_vol_signature)
      ray_vol_signature[idx] = RaySignature[]
  end
  
   
   return TriadTendencies(spec_grid, zeros(nxx, nyy, nzz, kpl, ml), ray_vol_signature,
      zeros(kpl, ml),
      kin_box, interp_coef)
end


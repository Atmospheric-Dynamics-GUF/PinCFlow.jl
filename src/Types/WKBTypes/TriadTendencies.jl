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
    St_k::B
    Col_int::B
    kin_box::KinematicBox
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

struct TriadTendencies{A <: AbstractArray{<: AbstractFloat, 5}}
    kp::AbstractVector{<:AbstractFloat}
    m::AbstractVector{<:AbstractFloat}
    wavespectrum::A
    st_k::A
    col_int::A
    kin_box::KinematicBox
end


function TriadTendencies(namelists::Namelists,
   domain::Domain,
   wkb_mode::NoWKB)::TriadTendencies

   println("Trying call kinematic box for NoWKB")
   
  kin_box = KinematicBox(wkb_mode)

  println("Succesfully called kinematic for NoWKB")
  
  return TriadTendencies(zeros(0), zeros(0), [zeros(0, 0, 0, 0, 0) for i in 1:3]...,
    kin_box)
end

function TriadTendencies(namelists::Namelists,
   domain::Domain,
   constants::Constants,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})::TriadTendencies

   println("Trying call kinematic box for WKB mode")

   (; lref) = constants
   (; kp_size, m_size) = namelists.triad

    # Non-dimensionalize domain boundaries.
    lx = namelists.domain.lx / lref
    ly = namelists.domain.ly / lref
    lz = namelists.domain.lz / lref
    lkp = namelists.triad.lkp * lref
    lm = namelists.triad.lm * lref

   (; nxx, nyy, nzz) = domain

   (; kp_size, m_size) = namelists.triad
  
   # Compute the grid in kp and m direction

   kpmin = sqrt((2 * pi / lx)^2 + (2 * pi / ly)^2)
   kpmax = lkp
   mmin = 2 * pi / lz
   mmax = lm
   kp = log_range(kpmin, kpmax, kp_size)
   m = log_range(mmin, mmax, m_size)
  

   println(kp)
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

   println("Succesfully called kinematic box for WKB mode")
   
   return TriadTendencies(kp, m,
      [zeros(nxx, nyy, nzz, kpl, ml) for i in 1:3]...,
      kin_box)
end


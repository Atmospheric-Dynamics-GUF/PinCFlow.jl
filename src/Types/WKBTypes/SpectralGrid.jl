struct SpectralGrid{A <: AbstractVector{<:AbstractFloat}, B <: Float64}
    kp::A
    m::A
    kpl::Integer
    ml::Integer
    lambdakp::B
    lambdam::B
    loglkp::B  #log of lambdakp
    loglm::B   #log of lambdam
end

function SpectralGrid(wkb_mode::NoWKB)::SpectralGrid
    
    return SpectralGrid(zeros(0, 0), zeros(0, 0), 0, 0, 0.0, 0.0, 0.0, 0.0)

end


function SpectralGrid(namelists::Namelists,
   constants::Constants,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})::SpectralGrid

   println("Trying call kinematic box for WKB mode")

   (; lref) = constants
   (; kp_size, m_size) = namelists.triad

    # Non-dimensionalize domain boundaries.
    lx = namelists.domain.lx / lref
    ly = namelists.domain.ly / lref
    lz = namelists.domain.lz / lref
    lkp = namelists.triad.lkp * lref
    lm = namelists.triad.lm * lref


   (; kp_size, m_size) = namelists.triad
  
   # Compute the grid in kp and m direction

   kpmin = sqrt((2 * pi / lx)^2 + (2 * pi / ly)^2)
   kpmax = lkp
   mmin = 2 * pi / lz
   mmax = lm
   kp = log_range(kpmin, kpmax, kp_size)
   m = log_range(mmin, mmax, m_size)
   kpl = length(kp)
   ml = length(m)

   return SpectralGrid(kp, m, kpl, ml, kp[2] / kp[1], m[2] / m[1], log(kp[2] / kp[1]), log(m[2] / m[1]))

    
end
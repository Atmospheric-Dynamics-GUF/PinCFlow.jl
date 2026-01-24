struct SpectralGrid{A <: AbstractVector{<:AbstractFloat}, B <: Float64}
    kp::A
    m::A
    kpc::A
    mc::A
    kpl::Integer
    ml::Integer
    lambdakp::B
    lambdam::B
    loglkp::B  #log of lambdakp
    loglm::B   #log of lambdam
end

function SpectralGrid(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad)::SpectralGrid
    
    return SpectralGrid(zeros(0), zeros(0), zeros(0), zeros(0), 0, 0, 0.0, 0.0, 0.0, 0.0)

end


function SpectralGrid(namelists::Namelists,
   constants::Constants,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
   triad_mode::Union{Triad2D, Triad3DIso})::SpectralGrid

   (; lref) = constants
   (; x_size, y_size, z_size) = namelists.domain
   (;k_size, l_size, m_size) = namelists.triad


    # Non-dimensionalize domain boundaries.
    lx = namelists.domain.lx / lref
    ly = namelists.domain.ly / lref
    lz = namelists.domain.lz / lref
    lk = namelists.triad.lk * lref
    ll = namelists.triad.ll * lref
    lm = namelists.triad.lm * lref

    #compute the grid in kp-direction

    if triad_mode == Triad2D() && k_size != 1 && x_size != 1
        kmin = 2 * pi / lx
        kmax = lk
        kp = log_range(kmin, kmax, k_size)
        lambdakp = kp[2] / kp[1]
        loglkp = log(lambdakp)
    elseif triad_mode == Triad2D() && l_size != 1 && y_size != 1
        lmin = 2 * pi / ly
        lmax = ll
        kp = log_range(lmin, lmax, l_size)
        lambdakp = kp[2] / kp[1]
        loglkp = log(lambdakp)
    elseif triad_mode == Triad3DIso() && k_size != 1 && y_size != 1 
        kpmin = sqrt((2 * pi / lx)^2 + (2 * pi / ly)^2)
        kpmax = sqrt(lk^2 + ll^2)
        kp = log_range(kpmin, kpmax, k_size)
        lambdakp = kp[2] / kp[1]
        loglkp = log(lambdakp)   
    else
        kp = zeros(1)
        lambdakp = 0.0
        loglkp = 0.0
    end


    #compute the grid in m-direction

    if m_size == 1 || z_size == 1
        m = zeros(1)
        lambdam = 0.0
        loglm = 0.0
        mc = zeros(0)
    else
        mmin = 2 * pi / lz
        mmax = lm
        m = log_range(mmin, mmax, m_size)
        lambdam = m[2] / m[1] 
        loglm = log(lambdam)
        mc = compute_edges(m)
        if triad_mode == Triad2D()   
            m = [-reverse(m); m]   #to include negative vertical wave number
            mc = [-reverse(mc); mc]
        end
    end
    kpc = compute_edges(kp)
    kpl = length(kp)
    ml = length(m)

   return SpectralGrid(kp, m, kpc, mc, kpl, ml, lambdakp, lambdam, loglkp, loglm)

    
end
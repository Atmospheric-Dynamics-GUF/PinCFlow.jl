struct SpectralGrid{A <: AbstractVector{<:AbstractFloat}, B <: Float64}
    kp::A  #grid center
    m::A
    kpc::A  #grid edges
    mc::A
    delkp::A
    delm::A
    kpl::Integer
    ml::Integer
    lambdakp::B
    lambdam::B
    loglkp::B  #log of lambdakp
    loglm::B   #log of lambdam
end

function SpectralGrid(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad)::SpectralGrid
    
    return SpectralGrid(zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), 0, 0, 0.0, 0.0, 0.0, 0.0)

end


function SpectralGrid(namelists::Namelists,
   constants::Constants,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
   triad_mode::Union{Triad2D, Triad3DIso})::SpectralGrid

   (; lref) = constants
   (; x_size, y_size, z_size) = namelists.domain
   (;k_size, l_size, m_size, k_max, l_max, m_max, k_min, l_min, m_min) = namelists.triad


    # Non-dimensionalize domain boundaries.
    lx = namelists.domain.lx / lref
    ly = namelists.domain.ly / lref
    lz = namelists.domain.lz / lref
    k_max = k_max * lref
    l_max = l_max * lref
    m_max = m_max * lref
    m_min = m_min * lref

    # Compute the grid in kp-direction.

    if triad_mode == Triad2D() && x_size == 1

        if wkb_mode != SingleColumn()
            error("Triad2D with x_size = 1 requires SingleColumn WKB mode.")
        end

        kmin = (k_min == 1) ? (2π / lx) : (k_min * lref)
        kmax = k_max

        dkp = 2π / lx

        nmin = ceil(Int, kmin / dkp - 100 * eps(Float64))
        nmax = floor(Int, kmax / dkp + 100 * eps(Float64))

        if nmin < 1 || nmax < nmin
            error("No discrete horizontal wave modes lie between k_min and k_max.")
        end

        kp = dkp .* collect(nmin:nmax)

        if length(kp) != k_size
            error("For x_size = 1, k_size must equal the number of discrete Fourier modes between k_min and k_max. Expected $(length(kp)), received $k_size.")
        end

        # k is discrete, so logarithmic-k quantities are not used.
        lambdakp = 0.0
        loglkp = 0.0

    elseif triad_mode == Triad2D() && k_size != 1 && x_size != 1

        kmin = (k_min == 1) ? (2π / lx) : (k_min * lref)
        kmax = k_max

        kp = log_range(kmin, kmax, k_size)
        lambdakp = kp[2] / kp[1]
        loglkp = log(lambdakp)

    elseif triad_mode == Triad2D() && l_size != 1 && y_size != 1

        lmin = (l_min == 1) ? (2π / lx) : (l_min * lref)
        lmax = l_max

        kp = log_range(lmin, lmax, l_size)
        lambdakp = kp[2] / kp[1]
        loglkp = log(lambdakp)

    elseif triad_mode == Triad3DIso() && k_size != 1 && y_size != 1

        kpmin = (k_min == 1 && l_min == 1) ?
            sqrt((2π / lx)^2 + (2π / ly)^2) :
            sqrt((k_min * lref)^2 + (l_min * lref)^2)

        kpmax = sqrt(k_max^2 + l_max^2)

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
        #mmin = 2 * pi / lz
        mmin = m_min
        mmax = m_max
        m = log_range(mmin, mmax, m_size)
        lambdam = m[2] / m[1] 
        loglm = log(lambdam)
        mc = compute_edges_centre(m)
        if triad_mode == Triad2D()   
            m = [-reverse(m); m]   #to include negative vertical wave number
            mc = [-reverse(mc); mc]
        end
    end
    kpc = compute_edges_centre(kp)
    kpl = length(kp)
    ml = length(m)
    delkp = zeros(kpl)
    delm = zeros(ml)

    if x_size == 1
    delkp .= 1.0
    else
        for kpi in eachindex(kp)
            delkp[kpi] = kpc[kpi + 1] - kpc[kpi]
        end
    end
    
    for mi in eachindex(m)
        if m[mi] > 0 
            delm[mi] = mc[mi + 2] - mc[mi + 1]
        else
            delm[mi] = mc[mi + 1] - mc[mi]
        end
    end
   return SpectralGrid(kp, m, kpc, mc, delkp, delm, kpl, ml, lambdakp, lambdam, loglkp, loglm)

    
end
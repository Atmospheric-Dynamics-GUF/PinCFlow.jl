"""
```julia
compute_spectral_cell_indices(
    state::State,
    xr::AbstractFloat,
    yr::AbstractFloat,
    dxr::AbstractFloat,
    dyr::AbstractFloat,
)::NTuple{4, <:Integer}
```

From the given spectral ray-volume position and extent, determine the indices of the spectral grid cells that contain the ray-volume edges and return them (in the order left, right, backward and forward).

# Arguments

  - `state`: Model state.

  - `kpr`: Ray-volume position in ``k_perp``.

  - `mrp`: Ray-volume position in ``m``.

  - `dkpr`: Ray-volume extent in ``k_perp``.

  - `dmr`: Ray-volume extent in ``m``.
"""
function compute_spectral_cell_indices end

function compute_spectral_cell_indices(
    state::State,
    kpr::AbstractFloat,
    mrr::AbstractFloat,
    dkpr::AbstractFloat,
    dmr::AbstractFloat,
)::NTuple{4, <:Integer}

    (; kp, m, kpc, mc, kpl, ml) = state.wkb.spec_tend.spec_grid
    mrp = abs(mrr)
    kp_l = kpr - dkpr/2
    kp_u = kpr + dkpr/2
    m_l = mrp - dmr/2
    m_u = mrp + dmr/2
    kpmin = kpmax = mmin = mmax = 0


    # Dealing for bounds if m includes negative wave numbers

    if m[1] < 0
        m0 = Int(ml / 2 + 1)
        m1 = ml
    else
        m0 = 1
        m1 = ml 
    end


    # Bounds check
    kp_lo = kpc[1]
    kp_hi = kpc[end]

    m_lo  = mc[m0]
    m_hi  = mc[end]

    out = false

    if kp_l < kp_lo
        println("Ray volume out of bounds in kp (lower)")
        println("  kp_l = ", kp_l, " < kp_min = ", kp_lo)
        out = true
    end

    if kp_u > kp_hi
        println("Ray volume out of bounds in kp (upper)")
        println("  kp_u = ", kp_u, " > kp_max = ", kp_hi)
        out = true
    end

    if m_l < m_lo
        if mrr >= 0
            println("Ray volume out of bounds in m (lower)")
            println("  m_l = ", m_l, " < m_min = ", m_lo)
            out = true
        else
            println("Ray volume out of bounds in -m (uper)")
            println("  -m_l = ", -m_l, " > -m_min = ", -m_lo)
            out = true
        end
    end

    if m_u > m_hi
        if mrr >= 0
            println("Ray volume out of bounds in m (upper)")
            println("  m_u = ", m_u, " > m_max = ", m_hi)
            out = true
        else
            println("Ray volume out of bounds in -m (lower)")
            println("  -m_u = ", -m_u, " < -m_max = ", -m_hi)
            out = true
        end
    end

    if out
        println("Ray-volume center and width:")
        println("  kpr = ", kpr, ", dkpr = ", dkpr)
        println("  mr  = ", mrr,  ", dmr  = ", dmr)
        error("Error: Ray volume out of spectral bound")
    end

    # ---- kp indices ----
    if kpl > 1
        for i in eachindex(kp)
            if kpc[i] <= kp_l
                kpmin = i
            end
            if kpc[i + 1] < kp_u
                kpmax = i + 1
            end
        end
    end

    # ---- m indices  ----
    if ml > 1
        for j in m0:m1
            if mc[j + 1] <= m_l
                mmin = j
            end
            if mc[j + 2] < m_u
                mmax = j + 1
            end
        end
    end

    kpmin = clamp(kpmin, 1, kpl)
    kpmax = clamp(kpmax, 1, kpl)
    mmin  = clamp(mmin,  m0, m1)
    mmax  = clamp(mmax,  m0, m1)

    #setting the correct indices for negative vertical wave number
    if mrr < 0
        (mmax, mmin) = (Int(ml - mmin + 1), Int(ml - mmax + 1))  
    end

    return (kpmin, kpmax, mmin, mmax)
end




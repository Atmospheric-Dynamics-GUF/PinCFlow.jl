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

  - `mr`: Ray-volume position in ``m``.

  - `dkpr`: Ray-volume extent in ``k_perp``.

  - `dmr`: Ray-volume extent in ``m``.
"""
function compute_spectral_cell_indices end

function compute_spectral_cell_indices(
    state::State,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    dkpr::AbstractFloat,
    dmr::AbstractFloat,
)::NTuple{4, <:Integer}
    (; kp_size, m_size) = state.namelists.triad
    (; kp, m) = state.wkb.spec_tend.spec_grid
    
    mr = abs(mr)
    kp_l = kpr - dkpr/2
    kp_u = kpr + dkpr/2
    m_l = mr - dmr/2
    m_u = mr + dmr/2
    kpmin = kpmax = mmin = mmax = 0

    # Bounds check
    """
    if (kp_l < kp[1] - (kp[2] - kp[1]) / 2) ||
       (kp_u > kp[end] + (kp[end] - kp[end-1]) / 2) ||
       (m_l < m[1] - (m[2] - m[1]) / 2) ||
       (m_u > m[end] + (m[end] - m[end-1]) / 2)
        error("Error: Ray volume out of spectral bound")
    end
    """
    kp_lo = kp[1]  - (kp[2]     - kp[1])     / 2
    kp_hi = kp[end] + (kp[end] - kp[end-1]) / 2

    m_lo  = m[1]   - (m[2]      - m[1])      / 2
    m_hi  = m[end] + (m[end]   - m[end-1])  / 2

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
        println("Ray volume out of bounds in m (lower)")
        println("  m_l = ", m_l, " < m_min = ", m_lo)
        out = true
    end

    if m_u > m_hi
        println("Ray volume out of bounds in m (upper)")
        println("  m_u = ", m_u, " > m_max = ", m_hi)
        out = true
    end

    if out
        println("Ray-volume center and width:")
        println("  kpr = ", kpr, ", dkpr = ", dkpr)
        println("  mr  = ", mr,  ", dmr  = ", dmr)
        error("Error: Ray volume out of spectral bound")
    end



    # ---- kp indices ----
    if kp_size > 1
        for i in eachindex(kp)
            if i == 1
                l_half = kp[i] - (kp[i+1] - kp[i]) / 2
                u_half = kp[i] + (kp[i+1] - kp[i]) / 2
            elseif i == kp_size
                l_half = kp[i] - (kp[i] - kp[i-1]) / 2
                u_half = kp[i] + (kp[i] - kp[i-1]) / 2
            else
                l_half = kp[i] - (kp[i] - kp[i-1]) / 2
                u_half = kp[i] + (kp[i+1] - kp[i]) / 2
            end
            if l_half <= kp_l
                kpmin = i
            end
            if u_half < kp_u
                kpmax = i + 1
            end
        end
    end

    # ---- m indices  ----
    if m_size > 1
        for j in eachindex(m)
            if j == 1
                l_half = m[j] - (m[j+1] - m[j]) / 2
                u_half = m[j] + (m[j+1] - m[j]) / 2
            elseif j == m_size
                l_half = m[j] - (m[j] - m[j-1]) / 2
                u_half = m[j] + (m[j] - m[j-1]) / 2
            else
                l_half = m[j] - (m[j] - m[j-1]) / 2
                u_half = m[j] + (m[j+1] - m[j]) / 2
            end
            if l_half <= m_l
                mmin = j
            end
            if u_half < m_u
                mmax = j + 1
            end
        end
    end

    return (kpmin, kpmax, mmin, mmax)
end


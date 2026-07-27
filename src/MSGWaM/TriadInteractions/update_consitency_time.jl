function update_consistency_time end

function update_consistency_time(
    consistency_time_min::Real,
    n_local::Real,
    dudz::Real,
    dndz::Real,
    kp::AbstractVector,
    m::AbstractVector,
    kpi_parent::Integer,
    mi_parent::Integer,
    kpi_1::Integer,
    mi_1::Integer,
    kpi_2::Integer,
    mi_2::Integer,
    eps_denom::Real,
    res_type::Sum
)::AbstractFloat

    if kpi_1 < 0 || kpi_2 < 0 || mi_1 < 0 || mi_2 < 0 
        return consistency_time_min
    end
    kp_parent = kp[kpi_parent]
    m_parent = m[mi_parent]
    kp_1 = kp[kpi_1]
    m_1 = m[mi_1]

    kp_2 = kp[kpi_2]
    m_2 = m[mi_2]

    abs(m_parent) < eps_denom && return consistency_time_min
    abs(m_1) < eps_denom && return consistency_time_min
    abs(m_2) < eps_denom && return consistency_time_min

    cz_parent = compute_cz(n_local, kp_parent, m_parent)
    cz_1 = compute_cz(n_local, kp_1, m_1)
    cz_2 = compute_cz(n_local, kp_2, m_2)

    
    d_delta_omega_dt =
        -(
            cz_1 * kp_1 +
            cz_2 * kp_2 -
            cz_parent * kp_parent
        ) * dudz -
        (
            cz_1 * kp_1 / abs(m_1) +
            cz_2 * kp_2 / abs(m_2) -
            cz_parent * kp_parent / abs(m_parent)
        ) * dndz

    if abs(d_delta_omega_dt) > eps_denom
        consistency_time_candidate =
            sqrt(2.0 * abs(1 / d_delta_omega_dt))

        if isfinite(consistency_time_candidate)
            consistency_time_min =
                min(consistency_time_min, consistency_time_candidate)
        end
    end

    return consistency_time_min
end

function update_consistency_time(
    consistency_time_min::Real,
    n_local::Real,
    dudz::Real,
    dndz::Real,
    kp::AbstractVector,
    m::AbstractVector,
    kpi_parent::Integer,
    mi_parent::Integer,
    kpi_1::Integer,
    mi_1::Integer,
    kpi_2::Integer,
    mi_2::Integer,
    eps_denom::Real,
    res_type::Difference
)::AbstractFloat

    if kpi_1 < 0 || kpi_2 < 0 || mi_1 < 0 || mi_2 < 0 
        return consistency_time_min
    end
    kp_parent = kp[kpi_parent]
    m_parent = m[mi_parent]
    kp_1 = kp[kpi_1]
    m_1 = m[mi_1]

    kp_2 = kp[kpi_2]
    m_2 = m[mi_2]

    abs(m_parent) < eps_denom && return consistency_time_min
    abs(m_1) < eps_denom && return consistency_time_min
    abs(m_2) < eps_denom && return consistency_time_min

    cz_parent = compute_cz(n_local, kp_parent, m_parent)
    cz_1 = compute_cz(n_local, kp_1, m_1)
    cz_2 = compute_cz(n_local, kp_2, m_2)

    
    d_delta_omega_dt =
        -(
            cz_1 * kp_1 -
            cz_2 * kp_2 -
            cz_parent * kp_parent
        ) * dudz -
        (
            cz_1 * kp_1 / abs(m_1) -
            cz_2 * kp_2 / abs(m_2) -
            cz_parent * kp_parent / abs(m_parent)
        ) * dndz

    if abs(d_delta_omega_dt) > eps_denom
        consistency_time_candidate =
            sqrt(2.0 * abs(1 / d_delta_omega_dt))

        if isfinite(consistency_time_candidate)
            consistency_time_min =
                min(consistency_time_min, consistency_time_candidate)
        end
    end

    return consistency_time_min
end
function update_dephasing_time end

function update_dephasing_time(
    dephasing_time_min::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    kp_parent::AbstractFloat,
    m_parent::AbstractFloat,
    kp_1::AbstractFloat,
    m_1::AbstractFloat,
    kp_2::AbstractFloat,
    m_2::AbstractFloat,
    eps_denom::AbstractFloat,
    ::Sum,
)::AbstractFloat

    abs(m_parent) < eps_denom && return dephasing_time_min
    abs(m_1) < eps_denom && return dephasing_time_min
    abs(m_2) < eps_denom && return dephasing_time_min

    cz_parent = compute_cz(n_local, kp_parent, m_parent)
    cz_1 = compute_cz(n_local, kp_1, m_1)
    cz_2 = compute_cz(n_local, kp_2, m_2)

    d_delta_omega_dt =
        -(cz_1 * kp_1 + cz_2 * kp_2 - cz_parent * kp_parent) * dudz -
        (cz_1 * kp_1 / abs(m_1) + cz_2 * kp_2 / abs(m_2) -
         cz_parent * kp_parent / abs(m_parent)) * dndz

    if isfinite(d_delta_omega_dt) && d_delta_omega_dt != 0.0
        dephasing_time_candidate = sqrt(2.0 / abs(d_delta_omega_dt))

        if isfinite(dephasing_time_candidate)
            dephasing_time_min = min(dephasing_time_min, dephasing_time_candidate)
        end
    end

    return dephasing_time_min
end

function update_dephasing_time(
    dephasing_time_min::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    kp_parent::AbstractFloat,
    m_parent::AbstractFloat,
    kp_1::AbstractFloat,
    m_1::AbstractFloat,
    kp_2::AbstractFloat,
    m_2::AbstractFloat,
    eps_denom::AbstractFloat,
    ::Difference,
)::AbstractFloat

    abs(m_parent) < eps_denom && return dephasing_time_min
    abs(m_1) < eps_denom && return dephasing_time_min
    abs(m_2) < eps_denom && return dephasing_time_min

    cz_parent = compute_cz(n_local, kp_parent, m_parent)
    cz_1 = compute_cz(n_local, kp_1, m_1)
    cz_2 = compute_cz(n_local, kp_2, m_2)

    d_delta_omega_dt =
        -(cz_1 * kp_1 - cz_2 * kp_2 - cz_parent * kp_parent) * dudz -
        (cz_1 * kp_1 / abs(m_1) - cz_2 * kp_2 / abs(m_2) -
         cz_parent * kp_parent / abs(m_parent)) * dndz

    if isfinite(d_delta_omega_dt) && d_delta_omega_dt != 0.0
        dephasing_time_candidate = sqrt(2.0 / abs(d_delta_omega_dt))

        if isfinite(dephasing_time_candidate)
            dephasing_time_min = min(dephasing_time_min, dephasing_time_candidate)
        end
    end

    return dephasing_time_min
end
function compute_st_k end

function compute_st_k(
    spec_tend::TriadTendencies,
    p::AbstractFloat,
    q::AbstractFloat,
    nk::AbstractFloat,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    compute_dephasing_time::Bool,
    triad_mode::Triad2D,
    res_type::Sum,
    )
    stk = 0.0
    dephasing_time_stk = Inf
    eps_denom = 1.0e-14

    kp1, kp2 = compute_kp1kp2(kpr, p, res_type)

    #----------------------------------------------------------
    # k = 1 + 2, branch +
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Sum())

    @ivy if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
            check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)

        if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
                i_m_2k1 = interaction_matrix(kp2, kpr, kp1, m2, mr, m1, Difference(), triad_mode)
                i_m_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kpr, mr, kp1, m1, kp2, m2, eps_denom, Sum())
                end
            end
        end
    end

    #----------------------------------------------------------
    # k = 1 + 2, branch -
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Difference())

    @ivy if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
            check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)

        if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
                i_m_2k1 = interaction_matrix(kp2, kpr, kp1, m2, mr, m1, Difference(), triad_mode)
                i_m_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kpr, mr, kp1, m1, kp2, m2, eps_denom, Sum())
                end
            end
        end
    end

    return stk, dephasing_time_stk
end

function compute_st_k(
    spec_tend::TriadTendencies,
    p::AbstractFloat,
    q::AbstractFloat,
    nk::AbstractFloat,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    compute_dephasing_time::Bool,
    triad_mode::Triad2D,
    res_type::Difference,
)
    stk = 0.0
    dephasing_time_stk = Inf
    eps_denom = 1.0e-14

    kp1, kp2 = compute_kp1kp2(kpr, q, res_type)

    #----------------------------------------------------------
    # 1 = k + 2, branch +
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Sum())

    @ivy if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
            check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)

        if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Sum(), triad_mode)
                i_m_21k = interaction_matrix(kp2, kp1, kpr, m2, m1, mr, Difference(), triad_mode)
                i_m_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k - n2 * n1 * i_m_k12) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kpr, mr, kp1, m1, kp2, m2, eps_denom, Difference())
                end
            end
        end
    end

    #----------------------------------------------------------
    # 1 = k + 2, branch -
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Difference())

    @ivy if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
            check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)

        if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Sum(), triad_mode)
                i_m_21k = interaction_matrix(kp2, kp1, kpr, m2, m1, mr, Difference(), triad_mode)
                i_m_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k - n2 * n1 * i_m_k12) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kpr, mr, kp1, m1, kp2, m2, eps_denom, Difference())
                end
            end
        end
    end

    # The two difference permutations are equivalent.
    stk *= 2.0

    return stk, dephasing_time_stk
end

# for x_size = 1
function compute_st_k(
    spec_tend::TriadTendencies,
    kp1i::Integer,
    kp2i::Integer,
    nk::AbstractFloat,
    kr::AbstractFloat,
    mr::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    compute_dephasing_time::Bool,
    triad_mode::Triad2D,
    res_type::Sum,
)
    (; kp) = spec_tend.spec_grid

    kp1 = kp[kp1i]
    kp2 = kp[kp2i]

    stk = 0.0
    dephasing_time_stk = Inf
    eps_denom = 1.0e-14

    #----------------------------------------------------------
    # k = 1 + 2, branch +
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kr, kp1, kp2, mr, Sum(), Sum())

    if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
       check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1i, m1, triad_mode)

        @ivy if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2i, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_k12 = interaction_matrix(kr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
                i_m_2k1 = interaction_matrix(kp2, kr, kp1, m2, mr, m1, Difference(), triad_mode)
                i_m_1k2 = interaction_matrix(kp1, kr, kp2, m1, mr, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kr, mr, kp1, m1, kp2, m2, eps_denom, Sum())
                end
            end
        end
    end

    #----------------------------------------------------------
    # k = 1 + 2, branch -
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kr, kp1, kp2, mr, Sum(), Difference())

    if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
       check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1i, m1, triad_mode)

        @ivy if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2i, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_k12 = interaction_matrix(kr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
                i_m_2k1 = interaction_matrix(kp2, kr, kp1, m2, mr, m1, Difference(), triad_mode)
                i_m_1k2 = interaction_matrix(kp1, kr, kp2, m1, mr, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kr, mr, kp1, m1, kp2, m2, eps_denom, Sum())
                end
            end
        end
    end

    return stk, dephasing_time_stk
end
# for x_size = 1
function compute_st_k(
    spec_tend::TriadTendencies,
    kp1i::Integer,
    kp2i::Integer,
    nk::AbstractFloat,
    kr::AbstractFloat,
    mr::AbstractFloat,
    n_local::AbstractFloat,
    dudz::AbstractFloat,
    dndz::AbstractFloat,
    compute_dephasing_time::Bool,
    triad_mode::Triad2D,
    res_type::Difference,
)
    (; kp) = spec_tend.spec_grid

    kp1 = kp[kp1i]
    kp2 = kp[kp2i]

    stk = 0.0
    dephasing_time_stk = Inf
    eps_denom = 1.0e-14

    #----------------------------------------------------------
    # 1 = k + 2, branch +
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kr, kp1, kp2, mr, Difference(), Sum())

    if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
       check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1i, m1, triad_mode)

        @ivy if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2i, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_1k2 = interaction_matrix(kp1, kr, kp2, m1, mr, m2, Sum(), triad_mode)
                i_m_21k = interaction_matrix(kp2, kp1, kr, m2, m1, mr, Difference(), triad_mode)
                i_m_k12 = interaction_matrix(kr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k - n2 * n1 * i_m_k12) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kr, mr, kp1, m1, kp2, m2, eps_denom, Difference())
                end
            end
        end
    end

    #----------------------------------------------------------
    # 1 = k + 2, branch -
    #----------------------------------------------------------

    m1, m2 = compute_m1m2(kr, kp1, kp2, mr, Difference(), Difference())

    if check_resolved_spectral_mode(spec_tend, kp1, m1, triad_mode) &&
       check_resolved_spectral_mode(spec_tend, kp2, m2, triad_mode)

        n1 = interpolate_nk(spec_tend, kp1i, m1, triad_mode)

        @ivy if nk != 0.0 || n1 != 0.0
            n2 = interpolate_nk(spec_tend, kp2i, m2, triad_mode)

            if (nk != 0.0 && n1 != 0.0) ||
               (nk != 0.0 && n2 != 0.0) ||
               (n1 != 0.0 && n2 != 0.0)

                i_p_1k2 = interaction_matrix(kp1, kr, kp2, m1, mr, m2, Sum(), triad_mode)
                i_m_21k = interaction_matrix(kp2, kp1, kr, m2, m1, mr, Difference(), triad_mode)
                i_m_k12 = interaction_matrix(kr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)

                dg = compute_g_prime(kp1, kp2, m1, m2)

                st_contribution = i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k - n2 * n1 * i_m_k12) / abs(dg)

                stk += st_contribution

                if compute_dephasing_time && st_contribution != 0.0
                    dephasing_time_stk = update_dephasing_time(dephasing_time_stk, n_local, dudz, dndz, kr, mr, kp1, m1, kp2, m2, eps_denom, Difference())
                end
            end
        end
    end

    # The two difference permutations are equivalent under
    # interchange of the partner modes.
    stk *= 2.0

    return stk, dephasing_time_stk
end

function compute_st_k(
    spec_tend::TriadTendencies,
    p::AbstractFloat,
    q::AbstractFloat,
    nk::AbstractFloat,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    triad_mode::Triad3DIso)::AbstractFloat

    # Integrand of the collision integral
    
    (kp1, kp2) = compute_kp1kp2(kpr, p, q)
                        
    # k = 1 + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Sum())
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1), triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, abs(m2), triad_mode) 
    vk12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk = vk12^2 * (n1 * n2 - nk * (n1 + n2)) / abs(dg)

    # k = 1 + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Difference())


    n1 = interpolate_nk(spec_tend, kp1, abs(m1), triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, abs(m2), triad_mode) 
    vk12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk += vk12^2 * (n1 * n2 - nk * (n1 + n2)) / abs(dg)

    # 1 = k + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Sum())

    
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1), triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, abs(m2), triad_mode) 
    v1k2 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode) 
    dg = compute_g_prime(kp1, m1, kp2, m2)
    stk -= 2 * v1k2^2 * (nk * n2 - n1 * (nk + n2)) / abs(dg)    

    # 1 = k + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Difference())

    
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1), triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, abs(m2), triad_mode) 
    v1k2 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode) 
    stk -= 2 * v1k2^2 * (nk * n2 - n1 * (nk + n2)) / abs(dg)    
    
    stk *= kp1 * kp2 
    
   

    return stk

end 
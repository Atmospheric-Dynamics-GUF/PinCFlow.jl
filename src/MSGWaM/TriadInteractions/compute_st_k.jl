function compute_st_k end

function compute_st_k(
    spec_tend::TriadTendencies,
    p::AbstractFloat,
    q::AbstractFloat,
    nk::AbstractFloat,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    nn::AbstractFloat,
    triad_mode::Triad2D)::AbstractFloat

    # Integrand of the collision integral
    
    (kp1, kp2) = compute_kp1kp2(kpr, p)
                        
    # k = 1 + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Sum())
    
    n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode) 
    i_p_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
    i_m_2k1 = interaction_matrix(kp2, kpr, kp1, m2, mr, m1, Difference(), triad_mode)
    i_m_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Difference(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk = i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)

    # k = 1 + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Difference())
    
    n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode) 
    i_p_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum(), triad_mode)
    i_m_2k1 = interaction_matrix(kp2, kpr, kp1, m2, mr, m1, Difference(), triad_mode)
    i_m_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Difference(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk += i_p_k12 * (n1 * n2 * i_p_k12 - nk * n1 * i_m_2k1 - nk * n2 * i_m_1k2) / abs(dg)


    # 1 = k + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Sum())
    
    n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode) 
    i_p_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Sum(), triad_mode)
    i_m_21k = interaction_matrix(kp2, kp1, kpr, m2, m1, mr, Difference(), triad_mode)
    i_m_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk -= 2 * i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k -  n2 * nk * i_m_k12) / abs(dg)    

    # 1 = k + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Difference())
    
    n1 = interpolate_nk(spec_tend, kp1, m1, triad_mode)
    n2 = interpolate_nk(spec_tend, kp2, m2, triad_mode) 
    i_p_1k2 = interaction_matrix(kp1, kpr, kp2, m1, mr, m2, Sum(), triad_mode)
    i_m_21k = interaction_matrix(kp2, kp1, kpr, m2, m1, mr, Difference(), triad_mode)
    i_m_k12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference(), triad_mode)
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk -= 2 * i_m_k12 * (nk * n2 * i_p_1k2 - n1 * nk * i_m_21k -  n2 * nk * i_m_k12) / abs(dg)   
    
    stk *= nn 
     
    return stk

end 

function compute_st_k end

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
function compute_st_k end

function compute_st_k(
    spec_tend::TriadTendencies,
    p::AbstractFloat,
    q::AbstractFloat,
    nk::AbstractFloat,
    kpr::AbstractFloat,
    mr::AbstractFloat)::AbstractFloat

    # Integrand of the collision integral
    
    (kp1, kp2) = compute_kp1kp2(kpr, p, q)
                        
    # k = 1 + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Sum())
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1))
    n2 = interpolate_nk(spec_tend, kp2, abs(m2)) 
    vk12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum())
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk = vk12^2 * (n1 * n2 - nk * (n1 + n2)) / abs(dg)

    # k = 1 + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Sum(), Difference())
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1))
    n2 = interpolate_nk(spec_tend, kp2, abs(m2)) 
    vk12 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Sum())
    dg = compute_g_prime(kp1, kp2, m1, m2)
    stk += vk12^2 * (n1 * n2 - nk * (n1 + n2)) / abs(dg)

    # 1 = k + 2, branch +
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Sum())
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1))
    n2 = interpolate_nk(spec_tend, kp2, abs(m2)) 
    v1k2 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference()) 
    dg = compute_g_prime(kp1, m1, kp2, m2)
    stk -= 2 * v1k2^2 * (nk * n2 - n1 * (nk + n2)) / abs(dg)    

    # 1 = k + 2, branch -
    (m1, m2) = compute_m1m2(kpr, kp1, kp2, mr, Difference(), Difference())
    
    n1 = interpolate_nk(spec_tend, kp1, abs(m1))
    n2 = interpolate_nk(spec_tend, kp2, abs(m2)) 
    v1k2 = interaction_matrix(kpr, kp1, kp2, mr, m1, m2, Difference()) 
    stk -= 2 * v1k2^2 * (nk * n2 - n1 * (nk + n2)) / abs(dg)    
    
    stk *= kp1 * kp2 
    
   

    return stk

end 
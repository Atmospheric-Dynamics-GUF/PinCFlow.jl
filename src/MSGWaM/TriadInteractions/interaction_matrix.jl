function interaction_matrix end

function interaction_matrix(
    kpr::AbstractFloat, 
    kp1::AbstractFloat,
    kp2::AbstractFloat,
    mr::AbstractFloat,
    m1::AbstractFloat,
    m2::AbstractFloat,
    res_type::Sum,
    triad_mode::Triad2D,
    )::AbstractFloat
 
<<<<<<< HEAD
    v_k12 = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *    #for the hydrostatic case
=======
    v_k12_p = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *    #for the hydrostatic case
>>>>>>> 9fc4c86c28f63b7940afe8f655ea9129a79106be
            (kp1 * m2 - kp2 * m1) * 
            (sign(kpr * kp2) * abs(mr * m2) - sign(kpr * kp1) * abs(mr * m1) - mr * m1 + mr * m2) /
            2 / sqrt(2)  / kpr / abs(m1 * m2)
    #omegar = compute_omega_hat(kpr, mr)
    #omega1 = compute_omega_hat(kp1, m1) 
    #omega2 = compute_omega_hat(kp2, m2) 
    #v_k12_p = - omega1 * omega2 * sqrt(omegar * omega1 * omega2) *
    #           (kp1 * m2 - kp2 * m1) *
    #           (kpr * kp2 / (omegar * omega2) - kpr * kp1 / (omegar * omega1) - mr * m1 + mr * m2) /
    #           2 / sqrt(2)  / kpr / kp1 / kp2    #it has to be noted that all apearance of N is avoided here and have been explicitly included in compute_st_k
<<<<<<< HEAD
    return v_k12
=======
    return v_k12_p
>>>>>>> 9fc4c86c28f63b7940afe8f655ea9129a79106be
end

function interaction_matrix(
    kpr::AbstractFloat, 
    kp1::AbstractFloat,
    kp2::AbstractFloat,
    mr::AbstractFloat,
    m1::AbstractFloat,
    m2::AbstractFloat,
    res_type::Difference,
    triad_mode::Triad2D,
    )::AbstractFloat
 
<<<<<<< HEAD
    v_k12 = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *
=======
    v_k12_m = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *
>>>>>>> 9fc4c86c28f63b7940afe8f655ea9129a79106be
            (kp1 * m2 - kp2 * m1) * 
            (sign(kpr * kp2) * abs(mr * m2) - sign(kpr * kp1) * abs(mr * m1) - mr * m1 - mr * m2) /
            2 / sqrt(2)  / kpr / abs(m1 * m2)  
    #omegar = compute_omega_hat(kpr, mr)
    #omega1 = compute_omega_hat(kp1, m1) 
    #omega2 = compute_omega_hat(kp2, m2) 
    #v_k12_m = - omega1 * omega2 * sqrt(omegar * omega1 * omega2) *
    #           (kp1 * m2 - kp2 * m1) *
    #           (kpr * kp2 / (omegar * omega2) - kpr * kp1 / (omegar * omega1) - mr * m1 - mr * m2 - 2) /
    #           2 / sqrt(2)  / kpr / kp1 / kp2 #it has to be noted that all apearance of N is avoided here and have been explicitly included in compute_st_k
<<<<<<< HEAD
    return v_k12         
=======
    return v_k12_m         
>>>>>>> 9fc4c86c28f63b7940afe8f655ea9129a79106be
end

function interaction_matrix(
    kpr::AbstractFloat, 
    kp1::AbstractFloat,
    kp2::AbstractFloat,
    mr::AbstractFloat,
    m1::AbstractFloat,
    m2::AbstractFloat,
    res_type::Sum,
    triad_mode::Triad3DIso,
    )::AbstractFloat

    kp1dotkp2 = (kpr^2 - kp1^2 - kp2^2)/2
    kprdotkp1 = (kpr^2 + kp1^2 - kp2^2)/2
    kprdotkp2 = (kpr^2 - kp1^2 + kp2^2)/2
    
    v_k12 = sqrt(kpr * kp1 * kp2 / 32) *
            (kp1dotkp2 * sqrt(abs(mr / m1 / m2)) / kp1 / kp2 +
            kprdotkp2 * sqrt(abs(m1 / mr / m2)) / kpr / kp2 +
            kprdotkp1 * sqrt(abs(m2 / mr / m1)) / kpr / kp1)

    return v_k12
end


function interaction_matrix(
    kpr::AbstractFloat, 
    kp1::AbstractFloat,
    kp2::AbstractFloat,
    mr::AbstractFloat,
    m1::AbstractFloat,
    m2::AbstractFloat,
    res_type::Difference,
    triad_mode::Triad3DIso,
    )::AbstractFloat

    kp1dotkp2 = (-kpr^2 + kp1^2 + kp2^2)/2
    kprdotkp1 = (kpr^2 + kp1^2 - kp2^2)/2
    kprdotkp2 = (-kpr^2 + kp1^2 - kp2^2)/2
    
    v_1k2 = sqrt(kpr * kp1 * kp2 / 32) *
            (kp1dotkp2 * sqrt(abs(mr / m1 / m2)) / kp1 / kp2 +
            kprdotkp2 * sqrt(abs(m1 / mr / m2)) / kpr / kp2 +
            kprdotkp1 * sqrt(abs(m2 / mr / m1)) / kpr / kp1)
            
    return v_1k2
end


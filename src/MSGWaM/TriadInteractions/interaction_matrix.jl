function interaction_matrix end



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
 
    v_k12 = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *
            (kp1 * m2 - kp2 * m1) * 
            (sign(kpr * kp2) * abs(mr * m2) - sign(kpr * kp1) * abs(mr * m1) - mr * m1 + mr * m2) /
            4  / kpr / abs(m1 * m2)            
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
    triad_mode::Triad2D,
    )::AbstractFloat
 
    v_k12 = - sqrt(kpr * kp1 * kp2 / (abs(mr * m1 * m2))) *
            (kp1 * m2 - kp2 * m1) * 
            (sign(kpr * kp2) * abs(mr * m2) - sign(kpr * kp1) * abs(mr * m1) - mr * m1 - mr * m2) /
            4  / kpr / abs(m1 * m2)            
    return v_k12
end


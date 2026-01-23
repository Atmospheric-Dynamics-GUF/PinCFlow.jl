function compute_g_prime end 

function compute_g_prime(
    kp1::AbstractFloat, 
    kp2::AbstractFloat,
    m1::AbstractFloat,
    m2::AbstractFloat
    )::AbstractFloat

    return (kp1 * sign(m1) / m1^2) - (kp2 * sign(m2) / m2^2)
 
end
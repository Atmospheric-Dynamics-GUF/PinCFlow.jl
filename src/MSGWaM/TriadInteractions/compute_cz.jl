function compute_cz end

function compute_cz(n_local::AbstractFloat, 
    kp_val::AbstractFloat, 
    m_val::AbstractFloat)::AbstractFloat

    return -n_local * kp_val * sign(m_val) / abs(m_val)^2
    
end
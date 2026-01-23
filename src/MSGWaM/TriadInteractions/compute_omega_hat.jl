function compute_omega_hat end

function compute_omega_hat(
    nn::AbstractFloat,
    kpr::AbstractFloat, 
    mr::AbstractFloat,
    )::AbstractFloat

    return nn * abs(kpr) / abs(mr)

end
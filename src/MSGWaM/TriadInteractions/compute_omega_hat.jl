function compute_omega_hat end

function compute_omega_hat(
    kpr::AbstractFloat, 
    mr::AbstractFloat,
    )::AbstractFloat

    return  abs(kpr) / sqrt(kpr^2 + mr^2)  #the multiplication of the bruint viasala frequency have been accounted in the prefactor in compute_st_k

end
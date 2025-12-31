function compute_kp1kp2 end

function compute_kp1kp2(
    kpr::AbstractFloat,
    p::AbstractFloat, 
    q::AbstractFloat,
    )::NTuple{2, <:AbstractFloat}

    kp1 = (kpr + p + q) / 2

    kp2 = (kpr - p + q) / 2

    return (kp1, kp2)

end
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

function compute_kp1kp2(
    kpr::AbstractFloat,
    p::AbstractFloat
    )::NTuple{2, <:AbstractFloat}

    kp1 = (kpr + p) / 2
    kp2 = (kpr - p) / 2
    if kp1 < 0 || kp2 < 0
        error("error in domain or kinematic box, kp1 or kp2 < 0")
    end
    return (kp1, kp2)

end
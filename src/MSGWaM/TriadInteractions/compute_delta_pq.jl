function compute_delta_pq end

function compute_delta_pq(
    kpr::AbstractFloat,
    p::AbstractFloat, 
    q::AbstractFloat,
    )::AbstractFloat

    return sqrt((kpr^2 - p^2) * (2 * kpr + q) * q ) / 2

end
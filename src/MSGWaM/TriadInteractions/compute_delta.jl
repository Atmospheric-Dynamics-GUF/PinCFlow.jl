function compute_delta end

function compute_delta(
    kpi::AbstractFloat,
    kpi1::AbstractFloat, 
    kpi2::AbstractFloat,
    )::AbstractFloat

    f1 = -kpi + kpi1 + kpi2
    f2 = kpi - kpi1 + kpi2
    f3 = kpi + kpi1 - kpi2
    f4 = kpi + kpi1 + kpi2

    return sqrt(f1 * f2 * f3 * f4) / 2

end
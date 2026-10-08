function compute_discrete_k_index end

function compute_discrete_k_index(
    kp::AbstractVector{<:AbstractFloat},
    kpvalue::AbstractFloat,
)::Union{Nothing, Integer}

    dkp = kp[2] - kp[1]

    # Outer boundaries of the represented discrete-k bins.
    kp_lo = kp[1] - dkp / 2
    kp_hi = kp[end] + dkp / 2

    tol = 1.0e-10 * max(abs(kp_hi), dkp)

    if kpvalue < kp_lo - tol || kpvalue > kp_hi + tol
        return nothing
    end

    # Assign the ray to the nearest discrete Fourier mode.
    kpi = round(Int, (kpvalue - kp[1]) / dkp) + 1

    return clamp(kpi, 1, length(kp))
end
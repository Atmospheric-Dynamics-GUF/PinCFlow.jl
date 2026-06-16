function compute_dphidz_center end

function compute_dphidz_center(
    phi::AbstractArray,
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    transform::Function,
)::AbstractFloat
    (; dz, jac) = state.grid

    jac_lower =
        2.0 * jac[i, j, k - 1] * jac[i, j, k] /
        (jac[i, j, k - 1] + jac[i, j, k])

    jac_upper =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])

    dz_lower = jac_lower * dz
    dz_upper = jac_upper * dz

    phi_lower = transform(phi[i, j, k - 1])
    phi_upper = transform(phi[i, j, k + 1])

    return (phi_upper - phi_lower) / (dz_lower + dz_upper)
end

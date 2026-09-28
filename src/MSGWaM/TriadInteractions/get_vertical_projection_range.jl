function get_vertical_projection_range end

function get_vertical_projection_range(
    state::State,
    iray::Integer,
    jray::Integer,
    k::Integer,
    zr::AbstractFloat,
    dzr::AbstractFloat,
)
    (; ko, k0) = state.domain

    if ko == 0 && k == k0 - 1
        return k0 - 1, k0 - 1
    end

    kmin = get_next_half_level(iray, jray, zr - dzr / 2, state; dkd = 1)
    kmax = get_next_half_level(iray, jray, zr + dzr / 2, state; dkd = 1)

    return kmin, kmax
end
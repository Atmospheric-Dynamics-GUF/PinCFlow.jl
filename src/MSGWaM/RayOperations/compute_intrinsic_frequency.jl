function compute_intrinsic_frequency(
    state::State,
    indices::NTuple{4, <:Integer},
)
    (; f_coriolis_dim) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; rays) = state.wkb

    zr = rays.z[indices...]
    kr = rays.k[indices...]
    lr = rays.l[indices...]
    mr = rays.m[indices...]
    khr = sqrt(kr^2 + lr^2)

    n2r = interpolate_stratification(zr, state, N2())
    fr = f_coriolis_dim * tref

    return branchr * sqrt(n2r * khr^2 + fr^2 * mr^2) / sqrt(khr^2 + mr^2)
end

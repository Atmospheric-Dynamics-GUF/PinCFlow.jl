function compute_intrinsic_frequency(
    state::State,
    indices::Ntuple{4, <:Integer},
)
    (; f_coriolis_dim) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; rays) = state.wkb

    zr = rays.z[indices...]
    wnrk = rays.k[indices...]
    wnrl = rays.l[indices...]
    wnrm = rays.m[indices...]
    wnrh = sqrt(wnrk^2 + wnrl^2)

    nnr = stratification(zr, 1)

    f_cor_nd = f_coriolis_dim * tref

    return branchr * sqrt(nnr * wnrh^2 + f_cor_nd^2 * wnrm^2) /
           sqrt(wnrh^2 + wnrm^2)
end

function get_chi_vertical_indices end

function get_chi_vertical_indices(state::State)
    (; i0, j0, k0, k1) = state.domain
    (; lref) = state.constants
    (; zc) = state.grid
    (; chi_z_min, chi_z_max) = state.namelists.triad

    z_min = min(chi_z_min, chi_z_max)
    z_max = max(chi_z_min, chi_z_max)

    # Single target height: use nearest model level.
    if chi_z_min == chi_z_max
        kk_target = k0
        distance_min = Inf

        for kk in k0:k1
            distance = abs(zc[i0, j0, kk] * lref - chi_z_min)

            if distance < distance_min
                distance_min = distance
                kk_target = kk
            end
        end

        return kk_target, kk_target
    end

    # Interval: use all cell centres inside it.
    kk_min = 0
    kk_max = 0

    for kk in k0:k1
        z = zc[i0, j0, kk] * lref

        if z >= z_min && z <= z_max
            if kk_min == 0
                kk_min = kk
            end

            kk_max = kk
        end
    end

    # Interval narrower than the grid spacing.
    if kk_min == 0
        z_target = 0.5 * (z_min + z_max)

        kk_target = k0
        distance_min = Inf

        for kk in k0:k1
            distance = abs(zc[i0, j0, kk] * lref - z_target)

            if distance < distance_min
                distance_min = distance
                kk_target = kk
            end
        end

        return kk_target, kk_target
    end

    return kk_min, kk_max
end
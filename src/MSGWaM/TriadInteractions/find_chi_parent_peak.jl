function find_chi_parent_peak end

function find_chi_parent_peak(
    was_pred,
    delkp,
    delm,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    kpi::Integer,
    mi_center::Integer,
    mi_radius::Integer,
    x_size::Integer,
)
    nm = size(was_pred, 5)

    mi_min = max(1, mi_center - mi_radius)
    mi_max = min(nm, mi_center + mi_radius)

    best_mi = 0
    best_action = 0.0

    @inbounds for mi in mi_min:mi_max
        action = get_chi_cell_action(
            was_pred,
            delkp,
            delm,
            ii,
            jj,
            kk,
            kpi,
            mi,
            x_size,
        )

        if action > best_action
            best_action = action
            best_mi = mi
        end
    end

    return best_mi, best_action
end
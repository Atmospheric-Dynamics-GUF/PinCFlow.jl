function track_chi_parent_modes! end

function track_chi_parent_modes!(
    state::State,
    ::Triad2D,
)
    (; i0, j0) = state.domain
    (; x_size, y_size) = state.namelists.domain
    (; output_variables) = state.namelists.output
    (; spec_tend) = state
    (; was_pred, chi_tracker) = spec_tend
    (; delkp, delm) = spec_tend.spec_grid

    if !(:chi_parent in output_variables)
        return nothing
    end

    if x_size != 1 || y_size != 1
        error(
            "chi_parent diagnostic currently requires x_size = 1 and y_size = 1.",
        )
    end

    ntracked = chi_tracker.ntracked[]

    if ntracked == 0
        return nothing
    end

    kk_min, kk_max = get_chi_vertical_indices(state)

    chi_tracker.kpi .= 0
    chi_tracker.mi .= 0
    chi_tracker.valid .= false

    anchor_kk = zeros(Int, ntracked)
    anchor_mi = zeros(Int, ntracked)

    # -------------------------------------------------------------------------
    # Identify an anchor for each component using its previous spectral peak.
    # -------------------------------------------------------------------------

    for p in 1:ntracked
        kpi = chi_tracker.kpi_ref[p]
        mi_ref = chi_tracker.mi_ref[p]

        if kpi == 0 || mi_ref == 0
            continue
        end

        mi_radius = max(1, chi_tracker.mi_radius[p])

        best_action = 0.0

        for kk in kk_min:kk_max
            mi, action = find_chi_parent_peak(
                was_pred,
                delkp,
                delm,
                i0,
                j0,
                kk,
                kpi,
                mi_ref,
                mi_radius,
                x_size,
            )

            if action > best_action
                best_action = action
                anchor_kk[p] = kk
                anchor_mi[p] = mi
            end
        end

        if anchor_kk[p] == 0
            continue
        end

        kk_anchor = anchor_kk[p]

        chi_tracker.kpi[p, kk_anchor] = kpi
        chi_tracker.mi[p, kk_anchor] = anchor_mi[p]
        chi_tracker.valid[p, kk_anchor] = true

        # ---------------------------------------------------------------------
        # Track upward in z by continuity.
        # ---------------------------------------------------------------------

        mi_previous = anchor_mi[p]

        for kk in (kk_anchor + 1):kk_max
            mi, action = find_chi_parent_peak(
                was_pred,
                delkp,
                delm,
                i0,
                j0,
                kk,
                kpi,
                mi_previous,
                mi_radius,
                x_size,
            )

            if action <= 0.0
                break
            end

            chi_tracker.kpi[p, kk] = kpi
            chi_tracker.mi[p, kk] = mi
            chi_tracker.valid[p, kk] = true

            mi_previous = mi
        end

        # ---------------------------------------------------------------------
        # Track downward in z by continuity.
        # ---------------------------------------------------------------------

        mi_previous = anchor_mi[p]

        for kk in (kk_anchor - 1):-1:kk_min
            mi, action = find_chi_parent_peak(
                was_pred,
                delkp,
                delm,
                i0,
                j0,
                kk,
                kpi,
                mi_previous,
                mi_radius,
                x_size,
            )

            if action <= 0.0
                break
            end

            chi_tracker.kpi[p, kk] = kpi
            chi_tracker.mi[p, kk] = mi
            chi_tracker.valid[p, kk] = true

            mi_previous = mi
        end
    end

    # -------------------------------------------------------------------------
    # If two initially distinct components are assigned to exactly the same
    # spectral cell at the same height, their identities are ambiguous there.
    # -------------------------------------------------------------------------

    for kk in kk_min:kk_max
        for p1 in 1:(ntracked - 1)
            if !chi_tracker.valid[p1, kk]
                continue
            end

            for p2 in (p1 + 1):ntracked
                if !chi_tracker.valid[p2, kk]
                    continue
                end

                same_cell =
                    chi_tracker.kpi[p1, kk] == chi_tracker.kpi[p2, kk] &&
                    chi_tracker.mi[p1, kk] == chi_tracker.mi[p2, kk]

                if same_cell
                    chi_tracker.valid[p1, kk] = false
                    chi_tracker.valid[p2, kk] = false
                end
            end
        end
    end

    # -------------------------------------------------------------------------
    # Update the persistent spectral reference for the next timestep.
    # -------------------------------------------------------------------------

    for p in 1:ntracked
        kk = anchor_kk[p]

        if kk == 0 || !chi_tracker.valid[p, kk]
            continue
        end

        chi_tracker.mi_ref[p] = chi_tracker.mi[p, kk]
    end

    return nothing
end

function track_chi_parent_modes!(
    state::State,
    ::Triad3DIso,
)
    return nothing
end
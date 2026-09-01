function diagnose_triad_timescales! end

function diagnose_triad_timescales!(
    state::State,
    ::Triad2D,
)
    (; i0, j0, master) = state.domain
    (; x_size, y_size) = state.namelists.domain
    (; output_variables) = state.namelists.output
    (; lref, tref) = state.constants
    (; zc) = state.grid
    (; n2) = state.atmosphere
    (; spec_tend) = state
    (; was_pred, col_int, chi_parent, chi_tracker) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; chi_action_rel_tol) = state.namelists.triad

    if !(:chi_parent in output_variables)
        return nothing
    end

    if x_size != 1 || y_size != 1
        error(
            "chi_parent diagnostic currently requires x_size = 1 and y_size = 1.",
        )
    end

    chi_parent .= NaN

    ntracked = chi_tracker.ntracked[]

    if ntracked == 0
        return nothing
    end

    if chi_action_rel_tol < 0.0 || chi_action_rel_tol >= 1.0
        error(
            "chi_action_rel_tol must satisfy 0 <= chi_action_rel_tol < 1.",
        )
    end

    kk_min, kk_max = get_chi_vertical_indices(state)

    # -------------------------------------------------------------------------
    # Find maximum contained action of each tracked component.
    # -------------------------------------------------------------------------

    action_max = zeros(Float64, ntracked)

    for p in 1:ntracked
        for kk in kk_min:kk_max
            if !chi_tracker.valid[p, kk]
                continue
            end

            kpi = chi_tracker.kpi[p, kk]
            mi = chi_tracker.mi[p, kk]

            action = get_chi_cell_action(
                was_pred,
                delkp,
                delm,
                i0,
                j0,
                kk,
                kpi,
                mi,
                x_size,
            )

            action_max[p] = max(action_max[p], action)
        end
    end

    chi_height = fill(NaN, ntracked)
    chi_kpi = zeros(Int, ntracked)
    chi_mi = zeros(Int, ntracked)

    # -------------------------------------------------------------------------
    # Maximum chi over the materially populated part of each tracked parent.
    # -------------------------------------------------------------------------

    for p in 1:ntracked
        if action_max[p] <= 0.0
            continue
        end

        action_threshold =
            chi_action_rel_tol * action_max[p]

        chi_max_p = -Inf

        for kk in kk_min:kk_max
            if !chi_tracker.valid[p, kk]
                continue
            end

            kpi = chi_tracker.kpi[p, kk]
            mi = chi_tracker.mi[p, kk]

            action = get_chi_cell_action(
                was_pred,
                delkp,
                delm,
                i0,
                j0,
                kk,
                kpi,
                mi,
                x_size,
            )

            if action <= action_threshold
                continue
            end

            was = was_pred[i0, j0, kk, kpi, mi]

            if was <= 0.0 || !isfinite(was)
                continue
            end

            omega =
                sqrt(n2[i0, j0, kk]) *
                compute_omega_hat(kp[kpi], m[mi])

            if omega == 0.0 || !isfinite(omega)
                continue
            end

            st = col_int[i0, j0, kk, kpi, mi]

            chi =
                2π * abs(st) /
                (abs(omega) * was)

            if !isfinite(chi)
                continue
            end

            if chi > chi_max_p
                chi_max_p = chi
                chi_height[p] = zc[i0, j0, kk] * lref
                chi_kpi[p] = kpi
                chi_mi[p] = mi
            end
        end

        if isfinite(chi_max_p)
            chi_parent[p] = chi_max_p
        end
    end

    # -------------------------------------------------------------------------
    # Diagnostic output.
    # -------------------------------------------------------------------------

    if master
        println("")
        println("WWT parent-mode diagnostic:")
        println("  Initialized wave modes        = ", length(chi_parent))
        println("  Tracked spectral components   = ", ntracked)

        if kk_min == kk_max
            println(
                "  Diagnostic height             = ",
                zc[i0, j0, kk_min] * lref / 1.0E3,
                " km",
            )
        else
            println(
                "  Diagnostic interval           = ",
                zc[i0, j0, kk_min] * lref / 1.0E3,
                " -- ",
                zc[i0, j0, kk_max] * lref / 1.0E3,
                " km",
            )
        end

        println("")

        for p in 1:ntracked
            println("  Tracked component ", p, ":")

            if isnan(chi_parent[p])
                println("    chi                          = NaN")
                println("")
                continue
            end

            println("    chi                          = ", chi_parent[p])
            println("    height at maximum chi        = ", chi_height[p] / 1.0E3, " km")
            println("    k                            = ", kp[chi_kpi[p]] / lref, " m^-1")
            println("    m                            = ", m[chi_mi[p]] / lref, " m^-1")
            println("    maximum contained action     = ", action_max[p])
            println("")
        end
    end

    return nothing
end

function diagnose_triad_timescales!(
    state::State,
    ::Triad3DIso,
)
    return nothing
end



#=function diagnose_triad_timescales end

function diagnose_triad_timescales(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    ::Triad2D;
    interaction_height::Real = 20.0E3,
    k_parents::Tuple{<:Real,<:Real} = (2π / 50.0E3, 2π / 50.0E3),
    m_signs::Tuple{<:Integer,<:Integer} = (+1, -1),
)
    (; master, k0, k1) = state.domain
    (; x_size) = state.namelists.domain
    (; lref, tref) = state.constants
    (; zc) = state.grid
    (; n2) = state.atmosphere
    (; spec_tend) = state
    (; wavespectrum, col_int, action_ref, chi_parent, chi_max) = spec_tend
    (; kp, m, delkp, delm) = spec_tend.spec_grid
    (; action_rel_tol) = state.namelists.triad

    #----------------------------------------------------------
    # Diagnose only at the physical level nearest the intended
    # interaction height.
    #----------------------------------------------------------

    kk_interaction = k0
    distance_min = Inf

    for kki in k0:k1
        distance = abs(zc[ii, jj, kki] * lref - interaction_height)

        if distance < distance_min
            distance_min = distance
            kk_interaction = kki
        end
    end

    if kk != kk_interaction
        return nothing
    end

    # Reset only when evaluating the selected interaction level.
    chi_parent .= NaN
    chi_max[] = NaN

    #----------------------------------------------------------
    # Locate each parent independently.
    #----------------------------------------------------------

    kpi_parent = zeros(Int, 2)
    mi_parent = zeros(Int, 2)
    was_parent_max = zeros(Float64, 2)

    for p in 1:2
        kp_target = k_parents[p] * lref

        if x_size == 1
            kpi_parent[p] = compute_discrete_k_index(kp, kp_target)
        else
            kpi_parent[p] = argmin(abs.(kp .- kp_target))
        end

        kpi = kpi_parent[p]

        for mi in eachindex(m)
            was = wavespectrum[ii, jj, kk, kpi, mi]

            if m_signs[p] * m[mi] > 0.0 && was > was_parent_max[p]
                was_parent_max[p] = was
                mi_parent[p] = mi
            end
        end
    end

    if any(iszero, mi_parent)
        if master
            println(
                "Timescale diagnostic: parent mode not found at z = ",
                zc[ii, jj, kk] * lref / 1.0E3,
                " km. mi_parent = ",
                mi_parent,
            )
        end

        return nothing
    end

    #----------------------------------------------------------
    # Check materially significant parent action.
    #----------------------------------------------------------

    action_floor = action_rel_tol * action_ref[]
    action_parent = zeros(Float64, 2)

    for p in 1:2
        kpi = kpi_parent[p]
        mi = mi_parent[p]

        spectral_measure =
            x_size == 1 ?
            abs(delm[mi]) :
            abs(delkp[kpi] * delm[mi])

        action_parent[p] =
            wavespectrum[ii, jj, kk, kpi, mi] *
            spectral_measure
    end

    if any(action_parent .<= action_floor)
        if master
            println(
                "Timescale diagnostic: parent action below threshold at z = ",
                zc[ii, jj, kk] * lref / 1.0E3,
                " km. action_parent = ",
                action_parent,
                ", action_floor = ",
                action_floor,
            )
        end

        return nothing
    end

    #----------------------------------------------------------
    # Parent-mode WWT parameters.
    #----------------------------------------------------------

    nn = sqrt(n2[ii, jj, kk])

    tau_lin = zeros(Float64, 2)
    tau_nl = zeros(Float64, 2)

    for p in 1:2
        kpi = kpi_parent[p]
        mi = mi_parent[p]

        was = wavespectrum[ii, jj, kk, kpi, mi]
        st = col_int[ii, jj, kk, kpi, mi]

        omega = nn * compute_omega_hat(kp[kpi], m[mi])

        tau_lin[p] =
            omega != 0.0 ?
            2π / abs(omega) :
            Inf

        tau_nl[p] =
            st != 0.0 ?
            abs(was / st) :
            Inf

        chi_parent[p] =
            isfinite(tau_nl[p]) && tau_nl[p] > 0.0 ?
            tau_lin[p] / tau_nl[p] :
            0.0
    end

    chi_max[] = max(chi_parent[1], chi_parent[2])

    #----------------------------------------------------------
    # Output.
    #----------------------------------------------------------

    if master
        println("")
        println("WWT parent-mode time-scale diagnostic:")
        println(
            "  Interaction level             = ",
            zc[ii, jj, kk] * lref / 1.0E3,
            " km",
        )
        println("")

        for p in 1:2
            kpi = kpi_parent[p]
            mi = mi_parent[p]

            was = wavespectrum[ii, jj, kk, kpi, mi]
            st = col_int[ii, jj, kk, kpi, mi]

            println("  Parent mode ", p, ":")
            println("    k                            = ", kp[kpi] / lref, " m^-1")
            println("    m_peak                       = ", m[mi] / lref, " m^-1")
            println("    wave-action density          = ", was)
            println("    collision integral           = ", st)
            println("    linear period                = ", tau_lin[p] * tref, " s")
            println("    linear period                = ", tau_lin[p] * tref / 3600.0, " h")
            println("    nonlinear time               = ", tau_nl[p] * tref, " s")
            println("    nonlinear time               = ", tau_nl[p] * tref / 3600.0, " h")
            println("    T_linear / T_nonlinear       = ", chi_parent[p])
            println("")
        end

        println("  Maximum parent-mode WWT parameter:")
        println("    chi_max                      = ", chi_max[])
        println("")
    end

    return nothing
end
=#
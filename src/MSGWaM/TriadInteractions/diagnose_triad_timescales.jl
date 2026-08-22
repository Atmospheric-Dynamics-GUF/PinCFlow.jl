function diagnose_triad_timescales end

function diagnose_triad_timescales(
    state::State,
    ii::Integer,
    jj::Integer,
    kk::Integer,
    ::Triad2D;
    interaction_height::AbstractFloat = 20.0E3,
    k_parent::AbstractFloat = 2π / 50.0E3,
)
    (; master, k0, k1) = state.domain
    (; x_size) = state.namelists.domain
    (; lref, tref) = state.constants
    (; zc) = state.grid
    (; n2) = state.atmosphere
    (; spec_tend) = state
    (; wavespectrum, col_int, action_ref) = spec_tend
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

    #----------------------------------------------------------
    # Horizontal spectral index of the parent waves.
    #
    # Spectral-grid wavenumbers are nondimensional, whereas
    # k_parent supplied above is dimensional [m^-1].
    #----------------------------------------------------------

    kp_parent = k_parent * lref

    if x_size == 1
        kpi = compute_discrete_k_index(kp, kp_parent)
    else
        kpi = argmin(abs.(kp .- kp_parent))
    end

    #----------------------------------------------------------
    # Locate the positive-m and negative-m parent peaks at the
    # parent horizontal wavenumber.
    #----------------------------------------------------------

    mi_positive = 0
    mi_negative = 0

    was_positive_max = 0.0
    was_negative_max = 0.0

    for mi in eachindex(m)
        was = wavespectrum[ii, jj, kk, kpi, mi]

        if m[mi] > 0.0 && was > was_positive_max
            was_positive_max = was
            mi_positive = mi
        elseif m[mi] < 0.0 && was > was_negative_max
            was_negative_max = was
            mi_negative = mi
        end
    end

    if mi_positive == 0 || mi_negative == 0
        return nothing
    end

    #----------------------------------------------------------
    # Make sure both parent peaks contain materially significant
    # action before interpreting their nonlinear time scales.
    #----------------------------------------------------------

    spectral_measure_1 =
        x_size == 1 ?
        abs(delm[mi_positive]) :
        abs(delkp[kpi] * delm[mi_positive])

    spectral_measure_2 =
        x_size == 1 ?
        abs(delm[mi_negative]) :
        abs(delkp[kpi] * delm[mi_negative])

    action_1 = was_positive_max * spectral_measure_1
    action_2 = was_negative_max * spectral_measure_2

    action_floor = action_rel_tol * action_ref[]

    if action_1 <= action_floor || action_2 <= action_floor
        return nothing
    end

    #----------------------------------------------------------
    # Local buoyancy frequency.
    #----------------------------------------------------------

    nn = sqrt(n2[ii, jj, kk])

    #----------------------------------------------------------
    # Parent mode 1: positive m
    #----------------------------------------------------------

    m1 = m[mi_positive]
    n1 = wavespectrum[ii, jj, kk, kpi, mi_positive]
    st1 = col_int[ii, jj, kk, kpi, mi_positive]

    omega1 = nn * compute_omega_hat(kp[kpi], m1)

    tau_lin_1 =
        omega1 != 0.0 ?
        2π / abs(omega1) :
        Inf

    tau_nl_1 =
        st1 != 0.0 ?
        abs(n1 / st1) :
        Inf

    #----------------------------------------------------------
    # Parent mode 2: negative m
    #----------------------------------------------------------

    m2 = m[mi_negative]
    n2_mode = wavespectrum[ii, jj, kk, kpi, mi_negative]
    st2 = col_int[ii, jj, kk, kpi, mi_negative]

    omega2 = nn * compute_omega_hat(kp[kpi], m2)

    tau_lin_2 =
        omega2 != 0.0 ?
        2π / abs(omega2) :
        Inf

    tau_nl_2 =
        st2 != 0.0 ?
        abs(n2_mode / st2) :
        Inf

    #----------------------------------------------------------
    # Convert nondimensional model times to physical seconds.
    #----------------------------------------------------------

    tau_lin_1_sec = tau_lin_1 * tref
    tau_lin_2_sec = tau_lin_2 * tref

    tau_nl_1_sec = tau_nl_1 * tref
    tau_nl_2_sec = tau_nl_2 * tref

    ratio_1 =
        isfinite(tau_nl_1) && tau_nl_1 > 0.0 ?
        tau_lin_1 / tau_nl_1 :
        0.0

    ratio_2 =
        isfinite(tau_nl_2) && tau_nl_2 > 0.0 ?
        tau_lin_2 / tau_nl_2 :
        0.0

    #----------------------------------------------------------
    # Output
    #----------------------------------------------------------

    if master
        println("")
        println("Intended ES parent-mode time-scale diagnostic:")
        println(
            "  Interaction level             = ",
            zc[ii, jj, kk] * lref / 1.0E3,
            " km",
        )
        println("")

        println("  Parent mode 1 (+m):")
        println("    k                            = ", kp[kpi] / lref, " m^-1")
        println("    m_peak                       = ", m1 / lref, " m^-1")
        println("    wave-action density          = ", n1)
        println("    collision integral           = ", st1)
        println("    linear period                = ", tau_lin_1_sec, " s")
        println("    linear period                = ", tau_lin_1_sec / 3600.0, " h")
        println("    nonlinear time               = ", tau_nl_1_sec, " s")
        println("    nonlinear time               = ", tau_nl_1_sec / 3600.0, " h")
        println("    T_linear / T_nonlinear       = ", ratio_1)
        println("")

        println("  Parent mode 2 (-m):")
        println("    k                            = ", kp[kpi] / lref, " m^-1")
        println("    m_peak                       = ", m2 / lref, " m^-1")
        println("    wave-action density          = ", n2_mode)
        println("    collision integral           = ", st2)
        println("    linear period                = ", tau_lin_2_sec, " s")
        println("    linear period                = ", tau_lin_2_sec / 3600.0, " h")
        println("    nonlinear time               = ", tau_nl_2_sec, " s")
        println("    nonlinear time               = ", tau_nl_2_sec / 3600.0, " h")
        println("    T_linear / T_nonlinear       = ", ratio_2)
        println("")
    end

    return nothing
end
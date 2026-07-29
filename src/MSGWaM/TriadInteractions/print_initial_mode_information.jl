function print_initial_mode_information end

function print_initial_mode_information(
    mode_information::AbstractVector,
    kp::AbstractVector,
    m::AbstractVector,
    lref::AbstractFloat,
    action_ref::AbstractFloat,
    support_tol::AbstractFloat,
    diagonal_connectivity::Bool,
)

    println("")
    println(repeat("-", 80))
    println("Initial spectral-mode identification")
    println(repeat("-", 80))
    println("Number of identified modes = ",
        length(mode_information))
    println("support_tol = ",
        support_tol)
    println("diagonal_connectivity = ",
        diagonal_connectivity)
    println("")

    for mode_index in eachindex(mode_information)

        mode =
            mode_information[mode_index]

        peak_kpi =
            mode.peak_index[1]

        peak_mi =
            mode.peak_index[2]

        println("Mode ", mode_index)
        println("  Number of occupied spectral cells = ",
            mode.cell_count)

        println(
            "  Peak spectral index = ",
            (peak_kpi, peak_mi),
        )

        println(
            "  Peak spectral coordinates = ",
            (
                kp[peak_kpi] / lref,
                m[peak_mi] / lref,
            ),
        )

        println(
            "  kpi index range = ",
            mode.kpi_min,
            ":",
            mode.kpi_max,
        )

        println(
            "  mi index range = ",
            mode.mi_min,
            ":",
            mode.mi_max,
        )

        println(
            "  kp range = [",
            kp[mode.kpi_min] / lref,
            ", ",
            kp[mode.kpi_max] / lref,
            "]",
        )

        println(
            "  m range = [",
            m[mode.mi_min] / lref,
            ", ",
            m[mode.mi_max] / lref,
            "]",
        )

        println(
            "  Peak contained action = ",
            mode.peak_action,
        )

        println("")
    end

    println(
        "action_ref = ",
        action_ref,
    )

    println(
        "action_ref corresponds to the weakest ",
        "identified initialized mode.",
    )

    println(repeat("-", 80))
    println("")

    return nothing
end
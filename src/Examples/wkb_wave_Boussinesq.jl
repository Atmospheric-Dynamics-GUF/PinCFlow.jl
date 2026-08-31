function wkb_wave_Boussinesq(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_wave_periodic.h5",
    plot_file::AbstractString = "wkb_wave_packed_periodic.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 1,
    y_size::Integer = 1,
    z_size::Integer = 100,
)
    lx = 10.0e3
    ly = 10.0e3
    lz = 10.0e3

    parameters = (
        k = 2 * pi / lz,
        l = 0.0,
        m = 20 * pi / lz,
        rx = 0.0,
        ry = 0.0,
        rz = 0.1,
        x0 = 0.0,
        y0 = 0.0,
        z0 = lz / 2,
        a0 = 0.1,
    )
    (; k, l, m) = parameters

    domain = DomainNamelist(;
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        x_size,
        y_size,
        z_size,
        vertical_boundary_condition = :Periodic,
    )

    atmosphere = AtmosphereNamelist(;
        model = :Boussinesq,
        background = :StableStratification,
    )

    state = State(Namelists(; atmosphere, domain))

    wkb = WKBNamelist(;
        wkb_mode = :MultiColumn,
        initial_wave_field = (alpha, x, y, z) -> (
            k,
            l,
            m,
            omega(state, parameters, x, y, z),
            wave_action_density(state, parameters, x, y, z),
        ),
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 60.0,
        output_variables = [:u],
        prepare_restart,
        tmax = 3600.0,
    )

    integrate(Namelists(; atmosphere, wkb, domain, output))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:u, 0.5, 0.5, 0.5, 2),
            display_figure,
            time_unit = :min,
        )
    end

    return
end

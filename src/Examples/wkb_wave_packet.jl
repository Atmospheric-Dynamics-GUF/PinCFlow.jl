# src/Examples/wkb_wave_packet.jl

function wkb_wave_packet(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_wave_packet.h5",
    plot_file::AbstractString = "wkb_wave_packet.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
    z0 = 10000

    lx = 20000
    ly = 20000
    lz = 30000

    parameters = (
        k = 8 * pi / lx,
        l = 8 * pi / ly,
        m = 8 * pi / lz,
        rx = 0.5,
        ry = 0.5,
        rz = 0.5,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 20000,
        a0 = 0.05,
    )
    (; k, l, m) = parameters

    model = :Compressible
    background = :LapseRates
    coriolis_frequency = 0.0001

    atmosphere = AtmosphereNamelist(; background, coriolis_frequency, model)

    domain = DomainNamelist(; lx, ly, lz, npx, npy, npz, x_size, y_size, z_size)

    grid = GridNamelist(; resolved_topography = (x, y) -> z0)

    output = OutputNamelist(;
        output_file,
        output_interval = 600,
        output_variables = [:dchidt0],
        prepare_restart,
        tmax = 600,
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

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        leading_order_impact = true,
        initial_chi = (x, y, z) -> z,
    )

    integrate(Namelists(; atmosphere, domain, grid, output, wkb, tracer))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:dchidt0, 0.5, 0.5, 0.5, 2);
            display_figure,
            time_unit = :min,
        )
    end

    return
end

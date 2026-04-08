# src/Examples/wkb_wave_packet.jl

function wkb_wave_packet(;
    x_size::Integer = 16,
    y_size::Integer = 16,
    z_size::Integer = 32,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_wave_packet.h5",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    plot_file::AbstractString = "examples/results/wkb_wave_packet.svg",
)
    lx = 20000.0
    ly = 20000.0
    lz = 40000.0

    parameters = (
        k = 16 * pi / lx,
        l = 16 * pi / ly,
        m = 32 * pi / lz,
        rx = 0.25,
        ry = 0.25,
        rz = 0.25,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 20000.0,
        a0 = 0.05,
    )
    (; k, l, m) = parameters

    model = :Compressible
    background = :Realistic
    coriolis_frequency = 0.0001

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    output = OutputNamelist(;
        output_file,
        save_ray_volumes = true,
        prepare_restart,
        output_interval = 900,
        tmax = 900,
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

    integrate(Namelists(; atmosphere, domain, output, wkb))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(plot_file, data, ("nr", 8, 8, 16, 2); time_unit = "min")
            return
        end
    end

    return
end

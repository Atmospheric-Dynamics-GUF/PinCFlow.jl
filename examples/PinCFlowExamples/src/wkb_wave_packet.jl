# examples/PinCFlowExamples/src/wkb_wave_packet.jl

function wkb_wave_packet(;
    x_size::Integer = 16,
    y_size::Integer = 16,
    z_size::Integer = 32,
    npx::Integer = 3,
    npy::Integer = 3,
    npz::Integer = 3,
    output_file::AbstractString = "wkb_wave_packet.h5",
    output_steps::Bool = false,
    save_ray_volumes::Bool = true,
    prepare_restart::Bool = false,
    visualize::Bool = true,
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

    model = Compressible()
    background = Realistic()
    coriolis_frequency = 0.0001

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    state = State(Namelists(; atmosphere, domain))

    output = OutputNamelist(;
        output_file,
        output_steps,
        output_interval = 900,
        tmax = 900,
        save_ray_volumes,
        prepare_restart,
    )

    wkb = WKBNamelist(;
        wkb_mode = MultiColumn(),
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
            plot_output(
                "examples/results/wkb_wave_packet.svg",
                data,
                ("nr", 8, 8, 16, 2);
                time_unit = "min",
            )
            return
        end
    end

    return
end

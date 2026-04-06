# PinCFlowExamples.jl/src/wkb_mountain_wave.jl

function wkb_mountain_wave(;
    x_size::Integer = 40,
    y_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::String = "wkb_mountain_wave.h5",
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
)
    h0 = 150.0
    l0 = 5000.0
    rl = 10
    rh = 2

    lx = 400000.0
    ly = 400000.0
    lz = 20000.0
    dxr = lx / 20
    dyr = ly / 20
    dzr = lz / 10
    alpharmax = 0.0179

    atmosphere = AtmosphereNamelist(;
        background = :LapseRates,
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    grid = GridNamelist(;
        resolved_topography = (x, y) ->
            x^2 + y^2 <= (rl * l0)^2 ?
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh /
            (rh + 1) : 0.0,
        unresolved_topography = (alpha, x, y) ->
            x^2 + y^2 <= (rl * l0)^2 ?
            (
                pi / l0,
                0.0,
                h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) / (rh + 1),
            ) : (0.0, 0.0, 0.0),
    )

    output = OutputNamelist(;
        output_file,
        save_ray_volumes = true,
        prepare_restart,
        output_steps,
    )

    sponge = SpongeNamelist(;
        lhs_sponge = (x, y, z, t, dt) ->
            alpharmax / 3 * (
                exp((abs(x) - lx / 2) / dxr) +
                exp((abs(y) - ly / 2) / dyr) +
                exp((z - lz) / dzr)
            ),
        relaxed_u = (x, y, z, t, dt) -> 10.0,
    )

    wkb = WKBNamelist(; wkb_mode = :MultiColumn)

    integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        h5open(output_file) do data
            plot_output(
                "PinCFlowExamples.jl/results/wkb_mountain_wave.svg",
                data,
                ("nr", 20, 20, 10, 2);
            )
            return
        end
    end

    return
end

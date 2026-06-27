# src/Examples/mountain_wave.jl

function mountain_wave(;
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "mountain_wave.h5",
    plot_file::AbstractString = "mountain_wave.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
)
    h0 = 100
    l0 = 1000

    lx = 20000
    ly = 20000
    lz = 5000

    dxr = lx / 2
    dyr = ly / 2
    dzr = lz / 2
    alpharmax = 0.0179

    atmosphere = AtmosphereNamelist(;
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(;
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        x_size = 20,
        y_size = 20,
        z_size = 20,
    )

    grid = GridNamelist(;
        resolved_topography = (x, y) -> h0 / (1 + (x^2 + y^2) / l0^2),
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 600,
        output_variables = [:w],
        prepare_restart,
        tmax = 600,
    )

    sponge = SpongeNamelist(;
        lhs_sponge = (x, y, z, t, dt) -> begin
            alpharx =
                abs(x) >= (lx - dxr) / 2 ?
                sin(pi * (abs(x) - (lx - dxr) / 2) / dxr)^2 : 0.0
            alphary =
                abs(y) >= (ly - dyr) / 2 ?
                sin(pi * (abs(y) - (ly - dyr) / 2) / dyr)^2 : 0.0
            alpharz =
                z >= lz - dzr ? sin(pi / 2 * (z - (lz - dzr)) / dzr)^2 : 0.0
            return alpharmax * (alpharx + alphary + alpharz) / 3
        end,
        relaxed_u = (x, y, z, t, dt) -> 10.0,
    )

    integrate(Namelists(; atmosphere, domain, grid, output, sponge))

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:w, 0.5, 0.5, 0.25, 2);
            time_unit = :min,
        )
    end

    return
end

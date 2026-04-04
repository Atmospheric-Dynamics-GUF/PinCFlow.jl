# PinCFlowExamples.jl/src/mountain_wave.jl

function mountain_wave(;
    x_size::Integer = 40,
    y_size::Integer = 40,
    z_size::Integer = 40,
    npx::Integer = 3,
    npy::Integer = 3,
    npz::Integer = 3,
    prepare_restart::Bool = false,
    output_steps::Bool = false,
    visualize::Bool = true,
)
    h0 = 100.0
    l0 = 1000.0

    lx = 20000.0
    ly = 20000.0
    lz = 20000.0
    dxr = lx / 2
    dyr = ly / 2
    dzr = lz / 2
    alpharmax = 0.0179

    atmosphere = AtmosphereNamelist(;
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    grid = GridNamelist(;
        resolved_topography = (x, y) -> h0 / (1 + (x^2 + y^2) / l0^2),
    )

    output = OutputNamelist(;
        output_file = "mountain_wave.h5",
        output_variables = (:w,),
        prepare_restart,
        output_steps,
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
        h5open("mountain_wave.h5") do data
            plot_output(
                "PinCFlowExamples.jl/results/mountain_wave.svg",
                data,
                ("w", 20, 20, 10, 2);
            )
            return
        end
    end

    return
end

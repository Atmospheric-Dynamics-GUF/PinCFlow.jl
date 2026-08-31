# src/Examples/wkb_mountain_wave.jl

"""
```julia
wkb_mountain_wave(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_mountain_wave.h5",
    plot_file::AbstractString = "wkb_mountain_wave.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
```

Run the WKB-mountain-wave example simulation.

# Keywords

The keywords are analogous to those of `cold_bubble`.
"""
function wkb_mountain_wave(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_mountain_wave.h5",
    plot_file::AbstractString = "wkb_mountain_wave.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 20,
    y_size::Integer = 20,
    z_size::Integer = 20,
)
    h0 = 300
    l0 = 5000
    rl = 10
    rh = 2

    lx = 200000
    ly = 200000
    lz = 5000

    dxr = lx / 20
    dyr = ly / 20
    dzr = lz / 10
    alpharmax = 0.0179

    atmosphere = AtmosphereNamelist(;
        coriolis_frequency = 0.0,
        initial_u = (x, y, z) -> 10.0,
    )

    domain = DomainNamelist(; lx, ly, lz, npx, npy, npz, x_size, y_size, z_size)

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
        output_interval = 300,
        output_variables = [:uw],
        prepare_restart,
        tmax = 300,
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

    turbulence = TurbulenceNamelist(; turbulence_scheme = :NoTurbulence)

    integrate(
        Namelists(; atmosphere, domain, grid, output, sponge, wkb, turbulence),
    )

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            plot_file,
            output_file,
            (:uw, 0.5, 0.5, 0.25, 2);
            display_figure,
            time_unit = :min,
        )
    end

    return
end

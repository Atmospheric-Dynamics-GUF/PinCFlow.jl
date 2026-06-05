# src/Examples/wkb_tke_wave_packet.jl

function wkb_tke_wave_packet(;
    x_size::Integer = 1,
    y_size::Integer = 1,
    z_size::Integer = 100,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_tke_wave_packet.h5",
)
    lx = 30e3
    ly = 30e3
    lz = 100e3

    parameters = (
        k = 2 * pi / lx,
        l = 2 * pi / ly,
        m = 2 * pi / 3e3,
        rx = 0.0,
        ry = 0.0,
        rz = 0.05,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 20e3,
        a0 = 1.0,
    )
    (; k, l, m) = parameters

    background = :Isothermal
    model = :Compressible
    coriolis_frequency = 1e-4

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    output = OutputNamelist(;
        output_file,
        output_variables = [:u, :v, :w, :rhop, :e],
        output_interval = 36,
        tmax = 3600,
    )

    state = State(Namelists(; atmosphere, domain))

    wkb = WKBNamelist(;
        wkb_mode = :MultiColumn,
        use_saturation = false,
        initial_wave_field = (alpha, x, y, z) -> (
            k,
            l,
            m,
            omega(state, parameters, x, y, z),
            wave_action_density(state, parameters, x, y, z),
        ),
        turbulent_damping = true,
    )

    integrate(Namelists(; atmosphere, domain, output, wkb))

    return
end

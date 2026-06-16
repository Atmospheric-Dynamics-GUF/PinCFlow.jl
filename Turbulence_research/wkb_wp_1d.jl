# Turbulence_research/wkb_wp_1d.jl

function wkb_wp_1d(;
    x_size::Integer = 1,
    y_size::Integer = 1,
    z_size::Integer = 100,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_wp_1d.h5",
    output_interval::AbstractFloat = 360.0,
    tmax::AbstractFloat = 3600.0,
    a0::AbstractFloat = 1.0,
    model::Symbol = :Compressible,
    background::Symbol = :Isothermal,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_setup::Symbol = :TracerOn,
)
    lx = 300.0e3
    ly = 30.0e3
    lz = 60.0e3

    parameters = (
        k = 2 * pi / lx,
        l = 0.0,
        m = 2 * pi / 3.0e3,
        rx = 0.0,
        ry = 0.0,
        rz = 0.05,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 30.e3,
        a0,
    )
    (; k, l, m) = parameters

    coriolis_frequency = 0.0

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    output = OutputNamelist(;
        output_file,
        save_ray_volumes = true,
        output_variables = [
            :u,
            :v,
            :w,
            :rhop,
            :tke,
            :shear_production,
            :buoyancy_production,
            :dtkedt,
            :gwshear,
            :gwbuoy,
            :e,
        ],
        output_interval,
        tmax,
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

    turbulence = TurbulenceNamelist(; momentum_coupling, entropy_coupling)

    tracer = TracerNamelist(; tracer_setup, initial_chi = (x, y, z) -> z)

    integrate(Namelists(; atmosphere, domain, output, wkb, turbulence, tracer))

    return
end

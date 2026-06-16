# Turbulence_research/wkb_wp_3d.jl

function wkb_wp(;
    x_size::Integer = 32,
    y_size::Integer = 1,
    z_size::Integer = 100,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    lx::AbstractFloat = 9000.e3,
    ly::AbstractFloat = 300.e3,
    lz::AbstractFloat = 100.e3,
    k::AbstractFloat = 0.0,
    l::AbstractFloat = 2 * pi / 300.e3,
    m::AbstractFloat = 2 * pi / 1.e3,
    rx::AbstractFloat = 1500.e3,
    ry::AbstractFloat = 0.0,
    rz::AbstractFloat = 5.e3,
    x0::AbstractFloat = 0.0,
    y0::AbstractFloat = 0.0,
    z0::AbstractFloat = 30.e3,
    a0::AbstractFloat = 1.0,
    version::AbstractFloat = 2.0,
    coriolis_frequency::AbstractFloat = 1.e-4,
    output_file::AbstractString = "wkb_wp_3d.h5",
    output_interval::AbstractFloat = 360.0,
    tmax::AbstractFloat = 3600.0,
    model::Symbol = :Compressible,
    background::Symbol = :Isothermal,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_setup::Symbol = :TracerOn,
    leading_order_impact::Bool = true,
    next_order_impact::Bool = true,
    turbulence_impact::Bool = true,
)
    parameters = (
        k = k,
        l = l,
        m = m,
        rx = rx,
        ry = ry,
        rz = rz,
        x0 = x0,
        y0 = y0,
        z0 = z0,
        a0 = a0,
        version = version,
    )

    atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

    domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)

    output = OutputNamelist(;
        output_file,
        output_variables = [
            :u,
            :v,
            :w,
            :rhop,
            :tke,
            :qchi,
            :dchidt0,
            :dchidtq,
            :dtkedt,
            :gwshear,
            :gwbuoy,
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

    tracer = TracerNamelist(;
        tracer_setup,
        initial_chi = (x, y, z) -> z,
        leading_order_impact,
        next_order_impact,
        turbulence_impact,
    )

    integrate(Namelists(; atmosphere, domain, output, wkb, turbulence, tracer))

    return
end

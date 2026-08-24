function wkb_wave_Boussinesq(;
    output_file::AbstractString = "wkb_wave_Boussinesq.h5",
    prepare_restart::Bool = false,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    x_size::Integer = 100,
    y_size::Integer = 1,
    z_size::Integer = 100,
    vertical_boundary_condition::Symbol = :Periodic,
    visualize::Bool = false,
    output_interval::AbstractFloat = 100.0,
    tmax::AbstractFloat = 100.0,
)
    lx = 400.0e3
    ly = 30.0e3
    lz = 10.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / ly,
        m = 20 * pi / lz,
        rx = 0.0,
        ry = 0.0,
        rz = 0.0,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 0.0,
        a0 = 0.7,
        version = 2,
    )
    (; k, l, m) = parameters

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 2.0e-3
    buoyancy_frequency = 2.0e-2

    domain = DomainNamelist(;
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        x_size,
        y_size,
        z_size,
        vertical_boundary_condition,
    )

    atmosphere = AtmosphereNamelist(;
        model,
        background,
        coriolis_frequency,
        buoyancy_frequency,
        kinematic_viscosity = 0.0,
        thermal_conductivity = 0.0,
        kinematic_diffusivity = 0.0,
    )

    state = State(Namelists(; atmosphere, domain))

    discretization = DiscretizationNamelist(; dtmax = 100)

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

    output = OutputNamelist(;
        output_file,
        output_interval,
        output_variables = [:rhop, :w, :v, :u, :tke, :dchidt0, :dchidtq],
        prepare_restart,
        tmax,
    )

    poisson = PoissonNamelist(; initial_cleaning = true)

    turbulence = TurbulenceNamelist(;
        turbulence_scheme = :TKEScheme,
        momentum_coupling = true,
        tracer_coupling = true,
    )

    kenv = 2 * pi / lx
    lenv = 0.0
    menv = 2 * pi / lz

    function chils(x::Real, y::Real, z::Real)
        return 1/2^3 *
               (1 + cos(kenv * x)) *
               (1 + cos(lenv * y)) *
               (1 + cos(menv * (z - lz / 2)))
    end

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) -> chils(x, y, z),
    )

    integrate(
        Namelists(;
            atmosphere,
            discretization,
            turbulence,
            wkb,
            domain,
            output,
            tracer,
            poisson,
        ),
    )

    return
end

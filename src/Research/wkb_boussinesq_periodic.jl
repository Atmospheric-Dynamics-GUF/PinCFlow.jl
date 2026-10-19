# src/Examples/wkb_boussinesq_periodic.jl

function wkb_boussinesq_periodic(;
    display_figure::Bool = true,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "wkb_boussinesq_periodic.h5",
    plot_file::AbstractString = "wkb_boussinesq_periodic.svg",
    prepare_restart::Bool = false,
    visualize::Bool = true,
    x_size::Integer = 40,
    y_size::Integer = 1,
    z_size::Integer = 20,
)
    lx = 400.0e3
    ly = 30.0e3
    lz = 20.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / 30.0e3,
        m = 2 * pi / 1.0e3,
        rx = 0.0,
        ry = 0.0,
        rz = 0.0,
        x0 = 0.0,
        y0 = 0.0,
        z0 = 0.0,
        a0 = 0.7,
        lz = lz,
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
        vertical_boundary_condition = :Periodic,
    )

    atmosphere = AtmosphereNamelist(;
        model,
        background,
        coriolis_frequency,
        buoyancy_frequency,
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
        turbulent_damping = true,
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 1.0,
        output_variables = [
            :dchidt0,
            :dchidtq,
            :tke,
            :u,
            :e,
            :dtkedt,
            :uu,
            :uv,
            :uw,
            :vv,
            :vw,
            :utheta,
            :vtheta,
            :uchi0,
            :vchi0,
            :wchi0,
            :qchi,
            :tracerdiffusion,
        ],
        prepare_restart,
        tmax = 1.0,
    )

    discretization = DiscretizationNamelist(; dtmax = 5.0)

    function chi_largescale(x::Real, z::Real)
        return 0.25 *
               (1 + cos(2 * pi * (z - lz / 2) / lz)) *
               (1 + cos(2 * pi * x / lx))
    end

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) -> chi_largescale(x, z),
        leading_order_impact = true,
        turbulence_wave_impact = true,
    )

    turbulence = TurbulenceNamelist(; gw_coupling = true)

    integrate(
        Namelists(;
            atmosphere,
            wkb,
            domain,
            output,
            turbulence,
            tracer,
            discretization,
        ),
    )

    if visualize && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        plot_output(
            "wkb_boussinesq_periodic-dchidt0.svg",
            output_file,
            (:dchidt0, 0.5, 0.5, 0.5, 61);
            display_figure,
            time_unit = :h,
        )
        plot_output(
            "wkb_boussinesq_periodic-dchidtq.svg",
            output_file,
            (:dchidtq, 0.5, 0.5, 0.5, 61);
            display_figure,
            time_unit = :h,
        )
    end

    return
end

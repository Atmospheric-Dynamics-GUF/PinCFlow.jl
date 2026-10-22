# src/Examples/boussinesq_periodic.jl

function boussinesq_periodic(;
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    output_file::AbstractString = "boussinesq_periodic.h5",
    x_size::Integer = 448,
    y_size::Integer = 32,
    z_size::Integer = 640,
    lambdah::AbstractFloat = 300.0e3,
    a0::AbstractFloat = 1.0,
)
    lx = 4000.0e3
    ly = 300.0e3
    lz = 20.0e3

    parameters = (
        k = 0.0,
        l = 2 * pi / lambdah,
        m = 2 * pi / 1.0e3,
        a0 = a0,
        lx = lx,
        ly = ly,
        lz = lz,
    )
    (; k, l, m) = parameters

    model = :Boussinesq
    background = :StableStratification
    coriolis_frequency = 1.0e-4
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
    (; g) = state.constants

    atmosphere = AtmosphereNamelist(;
        model,
        background,
        coriolis_frequency,
        buoyancy_frequency,
        buoyancy_initialization = :initial_thetap,
        initial_pip = (x, y, z) -> real(
            pihat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_thetap = (x, y, z) ->
            real(
                bhat(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ) / g * thetabar(state, x, y, z),
        initial_u = (x, y, z) -> real(
            uhat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_v = (x, y, z) -> real(
            vhat(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
        initial_w = (x, y, z) -> real(
            what(state, parameters, x, y, z) *
            exp(1im * phi(parameters, x, y, z)),
        ),
    )

    output = OutputNamelist(;
        output_file,
        output_interval = 1.0,
        output_variables = [:u, :v, :w, :tke, :rhop, :tracerdiffusion],
        tmax = 1.0,
    )

    function chi_largescale(z::Real)
        return (1 + cos(2 * pi * (z - lz / 2) / lz)) / 2
    end

    function dchidz(z::Real)
        return -pi / (lz) * sin(2 * pi * (z - lz / 2) / lz)
    end

    function chihat(
        state::State,
        parameters::NamedTuple,
        x::Real,
        y::Real,
        z::Real,
    )
        return -1im / omega(state, parameters, x, y, z) *
               (what(state, parameters, x, y, z) * dchidz(x, z))
    end

    tracer = TracerNamelist(;
        tracer_setup = :TracerOn,
        initial_chi = (x, y, z) ->
            chi_largescale(x, z) + real(
                chihat(state, parameters, x, y, z) *
                exp(1im * phi(parameters, x, y, z)),
            ),
    )

    lt = 30.0

    turbulence = TurbulenceNamelist(;
        initial_tke = (x, y, z) ->
            max(
                0,
                lt * (
                    lt * m^2 / 2 * (
                        abs(uhat(state, parameters, x, y, z))^2 +
                        abs(vhat(state, parameters, x, y, z))^2 - real(
                            (
                                uhat(state, parameters, x, y, z) *
                                uhat(state, parameters, x, y, z) +
                                vhat(state, parameters, x, y, z) *
                                vhat(state, parameters, x, y, z)
                            ) * exp(2im * phi(parameters, x, y, z)),
                        )
                    ) -
                    lt * (
                        n2(state, x, y, z) + real(
                            1im *
                            m *
                            bhat(state, parameters, x, y, z) *
                            exp(1im * phi(parameters, x, y, z)),
                        )
                    )
                ),
            ) / 2,
            momentum_coupling = false,
    )

    integrate(Namelists(; atmosphere, domain, output, tracer, turbulence))

    return
end

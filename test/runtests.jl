using JuliaFormatter
using PinCFlow
using Statistics
using Test
using HDF5
using MPI

Base.Experimental.@optlevel 1

configurations = Dict(
    Boussinesq() => Dict(
        UniformBoussinesq() => Dict(
            MountainWave() => "boussinesq_uniform_mountain_wave",
            WKBMountainWave() => "boussinesq_uniform_wkb_mountain_wave",
        ),
        StratifiedBoussinesq() => Dict(
            MountainWave() => "boussinesq_stratified_mountain_wave",
            WKBMountainWave() => "boussinesq_stratified_wkb_mountain_wave",
        ),
    ),
    PseudoIncompressible() => Dict(
        Isothermal() => Dict(
            MountainWave() => "pseudo_incompressible_isothermal_mountain_wave",
            WKBMountainWave() => "pseudo_incompressible_isothermal_wkb_mountain_wave",
            WavePacket() => "pseudo_incompressible_isothermal_wave_packet",
        ),
    ),
    Compressible() => Dict(
        Isothermal() => Dict(
            MountainWave() => "compressible_isothermal_mountain_wave",
            WKBMountainWave() => "compressible_isothermal_wkb_mountain_wave",
        ),
    ),
)

@testset "PinCFlow tests" begin
    for model in keys(configurations)
        for background in keys(configurations[model])
            for testcase in keys(configurations[model][background])
                title =
                    uppercasefirst(
                        replace(
                            configurations[model][background][testcase],
                            "_" => "-",
                            "wkb" => "WKB",
                        ),
                    ) * " tests"
                @testset "$title" begin
                    atmosphere = AtmosphereNamelist(;
                        specifyreynolds = false,
                        reinv = 0.0E+0,
                        mu_viscous_dim = 1.5E-5,
                        background = background,
                        buoyancy_frequency = 1.0E-2,
                        theta0_dim = 3.0E+2,
                        temp0_dim = 3.0E+2,
                        press0_dim = 1.0E+5,
                        backgroundflow_dim = (1.0E+1, 0.0E+0, 0.0E+0),
                        coriolis_frequency = 1.0E-4,
                    )

                    discretization = DiscretizationNamelist(;
                        cfl = 5.0E-1,
                        cfl_wave = 5.0E-1,
                        dtmin_dim = 1.0E-6,
                        dtmax_dim = 1.0E+3,
                        adaptive_time_step = true,
                        limitertype = MCVariant(),
                    )

                    domain = DomainNamelist(;
                        sizex = 5,
                        sizey = 5,
                        sizez = 5,
                        nbx = 3,
                        nby = 3,
                        nbz = 3,
                        lx_dim = 1.0E+5,
                        ly_dim = 1.0E+5,
                        lz_dim = 1.0E+4,
                        npx = 1,
                        npy = 1,
                        npz = 1,
                        base_comm = MPI.COMM_WORLD,
                    )

                    grid = GridNamelist(;
                        mountainheight_dim = 1.2E+3,
                        mountainwidth_dim = 5.0E+3,
                        mountain_case = 13,
                        height_factor = 2.0E+0,
                        width_factor = 1.0E+1,
                        spectral_modes = 1,
                        stretch_exponent = 1.2E+0,
                    )

                    ice = IceNamelist(; icesetup = NoIce())

                    output = OutputNamelist(;
                        output_variables = (),
                        save_ray_volumes = true,
                        prepare_restart = true,
                        restart = true,
                        iin = 1,
                        output_steps = false,
                        noutput = 1,
                        maxiter = 1,
                        outputtimediff = 3.6E+3,
                        maxtime = 3.6E+3,
                        input_file = "test/" *
                                     configurations[model][background][testcase] *
                                     ".h5",
                        output_file = "test/pincflow_output.h5",
                    )

                    poisson = PoissonNamelist(;
                        tolpoisson = 1.0E-8,
                        maxiterpoisson = 1000,
                        preconditioner = true,
                        dtau = 1.0E+0,
                        maxiteradi = 2,
                        initialcleaning = true,
                        relative_tolerance = false,
                    )

                    setting = SettingNamelist(;
                        model = model,
                        testcase = testcase,
                        zboundaries = SolidWallBoundaries(),
                    )

                    sponge = SpongeNamelist(;
                        spongelayer = true,
                        sponge_uv = false,
                        spongeheight = 1.0E-1,
                        spongealphaz_dim = 1.0E-2,
                        spongealphaz_fac = 1.0E+0,
                        unifiedsponge = true,
                        lateralsponge = true,
                        spongetype = ExponentialSponge(),
                        spongeorder = 1,
                        cosmosteps = 1,
                        relax_to_mean = false,
                        perturbation_period = 3.6E+3,
                        perturbation_amplitude = 1.0E-1,
                        relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
                    )

                    tracer = TracerNamelist(; tracersetup = NoTracer())

                    turbulence =
                        TurbulenceNamelist(; turbulencesetup = NoTurbulence())

                    wavepacket = WavePacketNamelist(;
                        wavepacketdim = 3,
                        lambdax_dim = 1.0E+4,
                        lambday_dim = 1.0E+4,
                        lambdaz_dim = 1.0E+3,
                        x0_dim = 0.0E+0,
                        y0_dim = 0.0E+0,
                        z0_dim = 5.0E+3,
                        sigmax_dim = 1.0E+5,
                        sigmay_dim = 1.0E+5,
                        sigmaz_dim = 1.0E+4,
                        a0 = 1.0E+0,
                        branch = -1,
                    )

                    wkb = WKBNamelist(;
                        xrmin_dim = -5.0E+4,
                        xrmax_dim = 5.0E+4,
                        yrmin_dim = -5.0E+4,
                        yrmax_dim = 5.0E+4,
                        zrmin_dim = 0.0E+0,
                        zrmax_dim = 1.0E+4,
                        nrxl = 1,
                        nryl = 1,
                        nrzl = 1,
                        nrk_init = 1,
                        nrl_init = 1,
                        nrm_init = 1,
                        nray_fac = 4,
                        fac_dk_init = 1.0E-1,
                        fac_dl_init = 1.0E-1,
                        fac_dm_init = 1.0E-1,
                        branchr = -1,
                        merge_mode = ConstantWaveAction(),
                        nsmth_wkb = 2,
                        lsmth_wkb = true,
                        sm_filter = Shapiro(),
                        zmin_wkb_dim = 0.0,
                        lsaturation = true,
                        alpha_sat = 1.0E+0,
                        wkb_mode = MultiColumn(),
                        blocking = true,
                        long_threshold = 2.5E-1,
                        drag_coefficient = 1.0E+0,
                        nwm = 1,
                    )

                    namelists = Namelists(;
                        domain,
                        output,
                        setting,
                        discretization,
                        poisson,
                        atmosphere,
                        grid,
                        sponge,
                        wkb,
                        tracer,
                        ice,
                        turbulence,
                        wavepacket,
                    )

                    integrate(namelists)

                    data = h5open("test/pincflow_output.h5")
                    reference = h5open(
                        "test/" *
                        configurations[model][background][testcase] *
                        ".h5",
                    )

                    for key in keys(reference)
                        @test all(isapprox.(data[key], reference[key]))
                    end

                    close(data)
                    close(reference)

                    rm("test/pincflow_output.h5")
                end
            end
        end
    end

    @testset "Code-format tests" begin
        for target in ("docs/make.jl", "examples/", "src/", "test/")
            @test format(target)
        end
    end
end

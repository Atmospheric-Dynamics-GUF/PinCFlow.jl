#using Revise

include("../src/PinCFlow.jl")

using .PinCFlow
using HDF5

domain = DomainNamelist(;
	x_size = 8,
	y_size = 1,
	z_size = 8,
	nbx = 3,
	nby = 3,
	nbz = 3,
	lx  = 4.0E+4,
	ly  = 1.0E+4,
	lz = 1.5E+4,
	npx = 1,
	npy = 1,
	npz = 1,
)

output = OutputNamelist(;
	output_variables = (:rhop, :pip, :w, :u, :thetap, :n2, :rhobar, :thetabar, :n, :qv, :q, :iaux1, :iaux2, :iaux3, :wwp, :epp, :thp),
	prepare_restart = false,
	restart = false,
	iin = -1,
	output_steps = false,
	noutput = 1,
	maxiter = 1,
	outputtimediff = 1.0, # 3.6E+1, #E+3
	maxtime = 30.0, #3.6E+1, #E+3
	input_file = "./test/pincflow_input.h5",
	output_file = "./test/pincflow_output.h5",
)

setting = SettingNamelist(;
	model = PseudoIncompressible(),
	testcase = WKBMultipleWavePackets(),
)

discretization = DiscretizationNamelist(;
	cfl = 5.0E-1,
	dtmin_dim = 1.0E-5,
	dtmax_dim = 6.0E+1,
	adaptive_time_step = true,
	limitertype = MCVariant(),
)

poisson = PoissonNamelist(;
	tolpoisson = 1.0E-8,
	maxiterpoisson = 1000,
	preconditioner = true,
	dtau = 4.0E+0,
	maxiteradi = 2,
	initialcleaning = false,
	relative_tolerance = false,
)

atmosphere = AtmosphereNamelist(;
	specifyreynolds = false,
	reinv = 0.0E+0,
	mu_viscous_dim = 0.0E+0,
	background = Isothermal(),
	temp0_dim = 3.0E+2,
	ground_pressure = 1.0E+5,
	initial_wind = (0.0E+0, 0.0E+0, 0.0E+0),
	coriolis_frequency = 1.0E-4,
)

grid = GridNamelist(;
	#mountainheight_dim = 1.0E+3,
	mountainheight_dim = 0.0E+3,
	mountainwidth_dim = 1.0E+3,
	mountain_case = 3,
	height_factor = 1.0E+0,
	width_factor = 1.0E+0,
	spectral_modes = 1,
	stretch_exponent = 1.0E+0,
)
sponge = SpongeNamelist(;
	spongelayer = false,
	alpharmax = 1.79E-2,
	betarmax = 0.0E+0,
	lateralsponge = false,
	spongetype = SinusoidalSponge(),
	relax_to_mean = false,
	relaxation_wind = (0.0E+0, 0.0E+0, 0.0E+0),
)
ice = IceNamelist(;
	icesetup = IceOn(),
	dt_ice = 1.0,
	nscx = 10,
	nscy = 1,
	nscz = 20,
	cloudcover = CloudCoverOn(),
	#    large_scale_ice = true,
	parameterized_nucleation = true,
	#    parameterized_nucleation = false,
	#     parameterized_sgs_q = true,
	parameterized_sgs_q = false,
)
wkb = WKBNamelist(;
	zrmin = -2.0E+4,
	xrmax_dim = 2.0E+4,
	zrmin = -5.0E+3,
	yrmax_dim = 5.0E+3,
	zrmin = 0.0E+0,
	zrmax_dim = 1.5E+4,
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
	branchr = 1,
	merge_mode = ConstantWaveAction(),
	nsmth_wkb = 2,
	#lsmth_wkb = true,
	lsmth_wkb = false,
	sm_filter = Shapiro(),
	zmin_wkb_dim = 0.0,
	lsaturation = false,
	#lsaturation = true,
	alpha_sat = 1.0E+0,
	wkb_mode = MultiColumn(),
	#blocking = true,
	blocking = false,
	long_threshold = 2.5E-1,
	drag_coefficient = 1.0E+0,
	nwm = 3,
)

multiwavepackets = MultiWavePacketNamelist(; random_wavepackets = true, nwm = 3)
# multiwavepackets = MultiWavePacketNamelist(;
#     nwm = 2,
#     wavepacketdim = [1, 1],
#     lambdax_dim = [1.0E+4, 1.0E+4],
#     lambday_dim = [0.0E+0   , 0.0E+0   ],
#     lambdaz_dim = [-2.0E+3, 2.0E+3],
#     x0_dim = [5.0E+3, 5.0E+3],
#     y0_dim = [5.0E+3, 5.0E+3],
#     z0_dim = [7.0E+3, 9.0E+3],
#     sigmax_dim = [0.0E+3, 0.0E+3],
#     sigmay_dim = [0.0E+3, 0.0E+3],
#     sigmaz_dim = [6.0E+3, 6.0E+3],
#     a0 = [0.12E+0, 0.12E+0],
# )

namelists = Namelists(;
	domain = domain,
	output = output,
	setting = setting,
	discretization = discretization,
	poisson = poisson,
	atmosphere = atmosphere,
	grid = grid,
	sponge = sponge,
	ice = ice,
	wkb = wkb,
	#    wavepacket=wavepacket,
	multiwavepackets = multiwavepackets,
)

integrate(namelists)

#include("../examples/visualization/fast_plot_tjl04.jl")

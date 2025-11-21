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
	lx = 4.0E+4,
	ly = 1.0E+4,
	lz = 1.5E+4,
	npx = 2,
	npy = 1,
	npz = 2,
)

output = OutputNamelist(;
	output_variables = (:rhop, :pip, :w, :u, :thetap, :n2, :rhobar, :thetabar, :n, :qv, :q, :iaux1, :iaux2, :iaux3, :wwp, :epp, :thp),
	prepare_restart = false,
	restart = false,
	iin = -1,
	output_steps = false,
	output_interval = 1.0,
	tmax = 30.0, #3.6E+1, #E+3
	output_file = "./test/pincflow_output.h5",
)

discretization = DiscretizationNamelist(;
	cfl_number = 5.0E-1,
	dtmin = 1.0E-5,
	dtmax = 6.0E+1,
	adaptive_time_step = true,
	limiter_type = MCVariant(),
)
poisson = PoissonNamelist(;
	tolerance = 1.0E-8,
	poisson_iterations = 1000,
	preconditioner = true,
	dtau = 4.0E+0,
	initial_cleaning = false,
	tolerance_is_relative = false,
)
atmosphere = AtmosphereNamelist(;
	model = PseudoIncompressible(),
	coriolis_frequency = 1.0E-4,
	kinematic_viscosity = 0.0E+0,
	thermal_conductivity = 0.0E+0,
)

ice = IceNamelist(;
	icesetup = IceOn(),
	ice_test_case = WKBMultipleWavePackets(),
	dt_ice = 1.0,
	nscx = 10,
	nscy = 1,
	nscz = 20,
	cloudcover = CloudCoverOn(),
	parameterized_nucleation = true,
	parameterized_sgs_q = false,
)
wkb = WKBNamelist(;
	nrx = 1,
	nry = 1,
	nrz = 1,
	nrk = 1,
	nrl = 1,
	nrm = 1,
	multiplication_factor = 4,
	dkr_factor = 1.0E-1,
	dlr_factor = 1.0E-1,
	dmr_factor = 1.0E-1,
	branch = 1,
	merge_mode = ConstantWaveAction(),
	filter_order = 2,
	smooth_tendencies = false,
	filter_type = Shapiro(),
	use_saturation = false,
	#lsaturation = true,
	saturation_threshold = 1.0E+0,
	#zmin_wkb_dim = 0.0,
	#lsaturation = false,
	#lsaturation = true,
	#alpha_sat = 1.0E+0,
	wkb_mode = MultiColumn(),
	#blocking = true,
	blocking = false,
	long_threshold = 2.5E-1,
	drag_coefficient = 1.0E+0,
	wave_modes = 3,
)

multiwavepackets = MultiWavePacketNamelist(; random_wavepackets = true, nwm = 3)

namelists = Namelists(;
	domain = domain,
	output = output,
	discretization = discretization,
	poisson = poisson,
	atmosphere = atmosphere,
	ice = ice,
	wkb = wkb,
	multiwavepackets = multiwavepackets,
)

integrate(namelists)

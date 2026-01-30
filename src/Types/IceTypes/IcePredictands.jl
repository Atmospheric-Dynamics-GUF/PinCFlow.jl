"""
```julia
IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for prognostic ice variables.

```julia
IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	variables::Variables,
)::IcePredictands
```

Construct an `IcePredictands` instance with dimensions and initial values depending on the general configuration of ice physics, by dispatching to the appropriate method.

```julia
IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::NoIce,
	variables::Variables,
)::IcePredictands
```

Construct an `IcePredictands` instance with zero-size arrays for configurations without ice physics.

```julia
IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::AbstractIce,
	variables::Variables,
)::IcePredictands
```

Construct an `IcePredictands` instance with all arrays initialized as ``z \\rho`` (non-dimensionalized).

# Fields

  - `n::A`: Ice-crystal number concentration.

  - `q::A`: Ice mixing ratio.

  - `qv::A`: Water-vapor mixing ratio.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `ice_setup`: General ice-physics configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.
"""
struct IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
	n::A
	q::A
	qv::A
end

function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	variables::Variables,
	iceconstants::IceConstants,
	iceforcing::IceForcing, # added
)
	(; ice_setup) = namelists.ice

	return IcePredictands(
		namelists,
		constants,
		domain,
		atmosphere,
		grid,
		ice_setup,
		variables,
		iceconstants,
		iceforcing, # added
	)
end

function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::NoIce,
	variables::Variables,
	iceconstants::IceConstants,
)
	n = zeros(0, 0, 0)
	q = zeros(0, 0, 0)
	qv = zeros(0, 0, 0)

	return IcePredictands(n, q, qv)
end

function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::AbstractIce,
	variables::Variables,
	iceconstants::IceConstants,
	iceforcing::IceForcing,
)
	(; nxx, nyy, nzz) = domain
	(; i0, i1, j0, j1, k0, k1) = domain
	(; x, zc) = grid
	(; rhobar, pbar) = atmosphere
	(; rho, rhop, u, v, w, pip, p) = variables.predictands
	(; kappainv, pref, gamma, lref) = constants
	(; ground_pressure) = namelists.atmosphere
	(; epsil0hat, meanMassIce, mRef) = iceconstants

	# **********************
	# ISSRegion as IC
	# **********************    

	# center ISSR
	z0_issr = 8.e3 # [m]
	# vertical width ISSR (standard deviation of gaussian dist.)
	sig_issr = 4.e3 # [m]
	# max value S in ISSR
	S_issr = 1.45 #iceconstants.S_c 

	# **********************

	n = zeros(nxx, nyy, nzz)
	q = zeros(nxx, nyy, nzz)
	qv = zeros(nxx, nyy, nzz)
	qv_profile = zeros(nzz)

	#nondim.
	z0_issr = z0_issr / lref
	sig_issr = sig_issr / lref

	#define upper/lower bounds of ISSR
	zMin_issr = z0_issr - sig_issr
	zMax_issr = z0_issr + sig_issr

	p0 = ground_pressure / pref
	S0_ini = 0

	qv_profile = zeros(nzz)

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
		#exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)
		exn_pMean = (pbar[i, j, k] / p0) ^ (gamma - 1)

		#pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		presMean = p0 * exn_pMean ^ kappainv #kappaInv = c_p/R

		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean

		#theta = pbar[i, j, k] / rho_full
		thetaMean = pbar[i, j, k] / rhoMean

		#temp = theta * exn_p
		tempMean = thetaMean * exn_pMean

		#psi = psat_ice(temp, iceconstants)
		psiMean = psat_ice(tempMean, iceconstants)

		#dimensional IC for n, q_v, q
		n0 = 0.0 #[kg**-1]
		#n0 = i

		#use S0 to define q_v
		#  !to be consistent with RayTracer simulation
		#  !NB: S0 is not initial S_i
		#  !S0  = q_v(0)* \mean p/ \mean p_si
		#  !S_i = q_v(0)* p/p_si

		#z_factor = 1.0 / 2.0 * (tanh( (zc[i, j, k] - zMin_issr) / (0.1 * sig_issr) ) - tanh( (zc[i, j, k] - zMax_issr) / (0.1 * sig_issr) ) )
		S0 = S_issr * exp(- (zc[i, j, k] - z0_issr) ^ 2 / 2.0 / sig_issr^2) # full gaussian profile

		#if ((zc[i, j, k] >= zMin_issr) && (zc[i, j, k] <= zMax_issr))
		#	S0 = S_issr * exp(- (zc[i, j, k] - z0_issr) ^ 2 / 2.0 / sig_issr^2)
		#else
		#	S0 = S0_ini
		#end

		qv0 = epsil0hat * S0 * psiMean / presMean # [kg/kg]

		q0 = meanMassIce * n0

		n[i, j, k] = rhoMean * n0 * mRef #\hat N = \hat \rho \hat n
		qv[i, j, k] = rhoMean * qv0
		q[i, j, k] = rhoMean * q0

		# set initial S to zero for testing
		#S0 = 0.0

		# gaussian profile of n for debugging centered at z0_issr and x=0
		#n[i, j, k] = rhoMean * mRef * 1.0e6 * exp(- (zc[i, j, k] - z0_issr) ^ 2 / 2.0 / sig_issr^2 / 4.0) * exp(- (x[i] ) ^ 2 / 2.0 / (sig_issr)^2 / 4.0) # [#/kg]
		#q[i, j, k] = n[i, j, k] * meanMassIce / mRef # [kg/kg]

	end

	# print max value in vertical profile of qv
	println("Max initial qv profile = ", maximum(qv[:,:,:]))
	println("Max initial n profile = ", maximum(n[:,:,:]))
	println("Max initial q profile = ", maximum(q[:,:,:]))

	qv_profile .= qv[i0, j0, :]

	# safe vertical profile of qv
	@views iceforcing.qv_ref .=  qv_profile

	# qv = 0 for debugging
	#qv .= 0.0

	#println("Initialized IcePredictands with ISSR as IC.")

	return IcePredictands(n, q, qv)
end
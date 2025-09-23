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
	icesetup::NoIce,
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
	icesetup::AbstractIce,
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

  - `icesetup`: General ice-physics configuration.

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
)
	(; icesetup) = namelists.ice

	return IcePredictands(
		namelists,
		constants,
		domain,
		atmosphere,
		grid,
		icesetup,
		variables,
		iceconstants,
	)
end

function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	icesetup::NoIce,
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
	icesetup::AbstractIce,
	variables::Variables,
	iceconstants::IceConstants,
)
	(; nxx, nyy, nzz) = domain
	(; i0, i1, j0, j1, k0, k1) = domain
	(; ztfc) = grid
	(; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = atmosphere
	(; rho, rhop, u, v, w, pip, p) = variables.predictands
	(; kappainv, pref, gamma, lref) = constants
	(; press0_dim) = namelists.atmosphere
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

	#nondim.
	z0_issr = z0_issr / lref
	sig_issr = sig_issr / lref

	#define upper/lower bounds of ISSR
	zMin_issr = z0_issr - sig_issr
	zMax_issr = z0_issr + sig_issr

	p0 = press0_dim / pref
	S0_ini = 0

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
		#exn_p = pip[i, j, k] + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
		exn_pMean = (pstrattfc[i, j, k] / p0) ^ (gamma - 1)

		#pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		presMean = p0 * exn_pMean ^ kappainv #kappaInv = c_p/R

		rhoMean = rhostrattfc[i, j, k]
		rho_full = rho[i, j, k] + rhoMean

		#theta = pstrattfc[i, j, k] / rho_full
		thetaMean = pstrattfc[i, j, k] / rhoMean

		#temp = theta * exn_p
		tempMean = thetaMean * exn_pMean

		#psi = psat_ice(temp, iceconstants)
		psiMean = psat_ice(tempMean, iceconstants)

		#dimensional IC for n, q_v, q
		n0 = 0.0 #[kg**-1]

		#use S0 to define q_v
		#  !to be consistent with RayTracer simulation
		#  !NB: S0 is not initial S_i
		#  !S0  = q_v(0)* \mean p/ \mean p_si
		#  !S_i = q_v(0)* p/p_si

		if ((ztfc[i, j, k] >= zMin_issr) && (ztfc[i, j, k] <= zMax_issr))
			S0 = S_issr * exp(- (ztfc[i, j, k] - z0_issr) ^ 2 / 2.0 / sig_issr^2)
		else
			S0 = S0_ini
		end

		qv0 = epsil0hat * S0 * psiMean / presMean # [kg/kg]

		q0 = meanMassIce * n0

		n[i, j, k] = rhoMean * n0 * mRef #\hat N = \hat \rho \hat n
		qv[i, j, k] = rhoMean * qv0
		q[i, j, k] = rhoMean * q0

	end

	return IcePredictands(n, q, qv)
end

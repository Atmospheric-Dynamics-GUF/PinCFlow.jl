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
	n_hom::A
	q_hom::A
	qv::A
	n_in::A # IN per m3
	n_het::A
	q_het :: A
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
	(; ice_setup, nucleation) = namelists.ice

	return IcePredictands(
		namelists,
		constants,
		domain,
		atmosphere,
		grid,
		ice_setup,
		nucleation,
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
	ice_setup::NoIce,
	nucleation::AbstractNucleation,
	variables::Variables,
	iceconstants::IceConstants,
)
	n_hom = zeros(0, 0, 0)
	q_hom = zeros(0, 0, 0)
	qv = zeros(0, 0, 0)
	n_in = zeros(0, 0, 0)
	n_het = zeros(0, 0, 0)
	q_het = zeros(0, 0, 0)

	return IcePredictands(n_hom, q_hom, qv, n_in, n_het, q_het)
end

# homogeneous nucleation
function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::AbstractIce,
	nucleation::HomogeneousOnly,
	variables::Variables,
	iceconstants::IceConstants,
)
	(; nxx, nyy, nzz) = domain
	
	n_hom = zeros(nxx, nyy, nzz)
	q_hom = zeros(nxx, nyy, nzz)
	qv = zeros(nxx, nyy, nzz)
	n_in = zeros(0, 0, 0)
	n_het = zeros(0, 0, 0)
	q_het = zeros(0, 0, 0)
	

	initialize_ice_fields!(n_hom, q_hom, qv, n_in, n_het, q_het, namelists, constants, domain, atmosphere, grid, ice_setup, nucleation, variables, iceconstants)

	return IcePredictands(n_hom, q_hom, qv, n_in, n_het, q_het)
end

# both nucleation
function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::AbstractIce,
	nucleation::BothNucleation,
	variables::Variables,
	iceconstants::IceConstants,
)
	(; nxx, nyy, nzz) = domain
	
	n_hom = zeros(nxx, nyy, nzz)
	q_hom = zeros(nxx, nyy, nzz)
	qv = zeros(nxx, nyy, nzz)
	n_in = zeros(nxx, nyy, nzz)
	n_het = zeros(nxx, nyy, nzz)
	q_het = zeros(nxx, nyy, nzz)

	initialize_ice_fields!(n_hom, q_hom, qv, n_in, n_het, q_het, namelists, constants, domain, atmosphere, grid, ice_setup, nucleation, variables, iceconstants)

	return IcePredictands(n_hom, q_hom, qv, n_in, n_het, q_het)
end

# heterogeneous nucleation
function IcePredictands(
	namelists::Namelists,
	constants::Constants,
	domain::Domain,
	atmosphere::Atmosphere,
	grid::Grid,
	ice_setup::AbstractIce,
	nucleation::HeterogeneousOnly,
	variables::Variables,
	iceconstants::IceConstants,
)
	(; nxx, nyy, nzz) = domain
	
	n_hom = zeros(0, 0, 0)
	q_hom = zeros(0, 0, 0)
	qv = zeros(nxx, nyy, nzz)
	n_in = zeros(nxx, nyy, nzz)
	n_het = zeros(nxx, nyy, nzz)
	q_het = zeros(nxx, nyy, nzz)
	

	initialize_ice_fields!(n_hom, q_hom, qv, n_in, n_het, q_het, namelists, constants, domain, atmosphere, grid, ice_setup, nucleation, variables, iceconstants)

	return IcePredictands(n_hom, q_hom, qv, n_in, n_het, q_het)
end


function initialize_ice_fields!(n_hom, q_hom, qv, n_in, n_het, q_het, namelists, constants, domain, atmosphere, grid, ice_setup, nucleation, variables, iceconstants)
	(; i0, i1, j0, j1, k0, k1) = domain
	(; x, zc) = grid
	(; rhobar, pbar) = atmosphere
	(; rho, rhop, u, v, w, pip, p) = variables.predictands
	(; kappainv, pref, gamma, lref) = constants
	(; ground_pressure) = namelists.atmosphere
	(; epsil0, meanMassIce, mRef, Nref) = iceconstants

	# **********************
	# ISSRegion as IC
	# **********************    

	# center ISSR
	z0_issr = 8.e3 # [m] # 8.e3
	# vertical width ISSR (standard deviation of gaussian dist.)
	sig_issr = 4.e3 # [m]
	# max value S in ISSR
	S_issr = 1.45 #iceconstants.S_c 
	#S_issr = 0.9 # het. nucleation

	#nondim.
	z0_issr = z0_issr / lref
	sig_issr = sig_issr / lref

	#define upper/lower bounds of ISSR
	zMin_issr = z0_issr - 1 * sig_issr
	zMax_issr = z0_issr + 1 * sig_issr

	# define horizontal bounds of ISSR (just for advection tests)
	xMax_issr = 5.e3 / lref
	xMin_issr = -5.e3 / lref

	p0 = ground_pressure / pref
	S0_ini = 0

	# advection test n_in
	x0 = 0.2 * maximum(x)
	sigma_x = 0.1 * maximum(x)
	z0 = z0_issr
	sigma_z = sig_issr
	A  = 1e5

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

		psiMean = psat_ice(tempMean, iceconstants)		
		
		if ((zc[i, j, k] >= zMin_issr) && (zc[i, j, k] <= zMax_issr))
			S0 = S_issr * exp(- (zc[i, j, k] - z0_issr) ^ 2 / 2.0 / sig_issr^2)
		else
			S0 = S0_ini
		end
		
		# --- prognostic ice variables ---
		n0 = 0.0 #[kg^-1]
		# non-dimensional IC for n, q_v, q
		qv0 = epsil0 * S0 * psiMean / presMean # [kg/kg] 
		q0 = meanMassIce * n0

		qv[i, j, k] = rhoMean * qv0

		if nucleation isa HomogeneousOnly || nucleation isa BothNucleation
			n_hom[i, j, k] = rhoMean * n0 * mRef
			q_hom[i, j, k] = rhoMean * q0
		end

		if nucleation isa HeterogeneousOnly || nucleation isa BothNucleation
			#n_in0 = 5.e7 # [kg^-1]
			#n_in[i,j,k] = n_in0 * mRef * rhoMean
			# Use a initial condition in units of 1/m^3 instead of 1/kg 
			n_in0 = 10.e3 # [1/m^3]
			n_in[i,j,k] = n_in0 * mRef
			n_het[i, j, k] = rhoMean * n0 * mRef
			q_het[i, j, k] = rhoMean * q0
		end
	end
end
"""
```julia
activate_orographic_source!(
	state::State,
	omi_ini::AbstractArray{<:AbstractFloat, 4},
	wnk_ini::AbstractArray{<:AbstractFloat, 4},
	wnl_ini::AbstractArray{<:AbstractFloat, 4},
	wnm_ini::AbstractArray{<:AbstractFloat, 4},
	wad_ini::AbstractArray{<:AbstractFloat, 4},
)
```

Compute ray-volume properties in the launch layer (i.e. at `k = k0 - 1`) for the initialization of MSGWaM.

Sets the launch-layer values of arrays for initial ray-volume properties (intrinsic frequencies, wavenumbers and wave-action densities). For this purpose, the horizontal components of the resolved wind, the background density and the squared buoyancy frequency are vertically averaged between the surface and an approximation for the summits of the unresolved orography. The vertical averages are then used to compute a non-dimensionalized mountain wave amplitude, from which an approximate reduction of the generated wave amplitude due to blocking is inferred (see below). Afterwards, the ray-volume properties are obtained by calling `compute_orographic_mode` with the correspondingly scaled mode of the orographic spectrum and the vertical averages as arguments.

```julia
activate_orographic_source!(state::State)
```

Launch ray volumes that represent unresolved orographic gravity waves.

In each column of MPI processes at the lower boundary, this method first computes vertical averages of the horizontal components of the resolved wind, the background density and the squared buoyancy frequency between ``h_\\mathrm{b}`` (the surface) and ``h_\\mathrm{b} + \\Delta h`` (an approximation for the summits of the unresolved orography, with ``\\Delta h = \\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right|``). The vertical averages are then used to compute a non-dimensionalized mountain wave amplitude, from which an approximate reduction of the generated wave amplitude due to blocking, as well as the depth of the blocked layer, is inferred. A loop over the spectral modes of the unresolved orography follows, in which the properties of each mode are computed, using `compute_orographic_mode` with the scaled mode of the orographic spectrum and vertical averages as arguments, and corresponding ray volumes are launched at `k = k0 - 1`.

The parameterization of blocking is built around the non-dimensionalized mountain wave amplitude, or Long number,

```math
\\mathrm{Lo} = \\frac{N_h \\Delta h}{\\left|\\boldsymbol{u}_h\\right|},
```

where ``N_h`` is the square root of the vertically averaged squared buoyancy frequency and ``\\boldsymbol{u}_h`` is the vertically averaged resolved horizontal wind. This number is used to estimate the depth of the blocked layer as

```math
\\Delta z_\\mathrm{B} = 2 \\Delta h \\max \\left(0, \\frac{\\mathrm{Lo} - C}{\\mathrm{Lo}}\\right),
```

where ``C`` is a critical value represented by the model parameter `state.namelists.wkb.long_threshold`. The corresponding scaling of the orographic spectrum is given by

```math
r \\left(\\mathrm{Lo}\\right) = \\frac{2 \\Delta h - \\Delta z_\\mathrm{B}}{2 \\Delta h} = \\min \\left(1, \\frac{C}{\\mathrm{Lo}}\\right),
```

so that ``\\Delta z_\\mathrm{B} = 2 \\Delta h \\left(1 - r\\right)``. In addition to the reduction of the mountain-wave amplitude, the present blocked-layer scheme adds a blocked-flow drag to the mean-flow impact. This is implemented in [`PinCFlow.MSGWaM.MeanFlowEffect.apply_blocked_layer_scheme!`](@ref).

The launch algorithm distinguishes between the following situations (regarding previously launched ray volumes).

 1. There is no ray volume with nonzero wave-action density. A new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density that has at least partially passed through the lower boundary. The ray volume is either clipped or extended, such that its lower edge coincides with the surface, and the part below the surface is discarded. Then, it is assigned to the first model layer `k0`, i.e. its indices are changed from `(iray, ix, jy, k0 - 1)` to `(jray, ix, jy, k0)`, where `jray` is the new last ray-volume index at `(ix, jy, k0)`. Finally, a new ray volume is launched.

 1. There is a ray volume with nonzero wave-action density, which has not yet crossed the lower boundary. It is replaced with a new one.

# Arguments

  - `state`: Model state.

  - `omi_ini`: Array for intrinsic frequencies.

  - `wnk_ini`: Array for zonal wavenumbers.

  - `wnl_ini`: Array for meridional wavenumbers.

  - `wnm_ini`: Array for vertical wavenumbers.

  - `wad_ini`: Array for wave-action densities.

# See also

  - [`PinCFlow.MSGWaM.RaySources.compute_orographic_mode`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function activate_multiplewavepackets_source! end

function activate_multiplewavepackets_source!(
	state::State,
	omi_ini::AbstractArray{<:AbstractFloat, 4},
	wnk_ini::AbstractArray{<:AbstractFloat, 4},
	wnl_ini::AbstractArray{<:AbstractFloat, 4},
	wnm_ini::AbstractArray{<:AbstractFloat, 4},
	wad_ini::AbstractArray{<:AbstractFloat, 4},
)
	(; coriolis_frequency) = state.namelists.atmosphere
	(; branchr) = state.namelists.wkb
	(; tref, lref) = state.constants
	(; io, jo, ko, i0, i1, j0, j1, k0, k1) = state.domain
	(; ztfc, x, y) = state.grid
	(; random_wavepackets, wavepacketdim, lambdax_dim, lambday_dim, lambdaz_dim,
		x0_dim, y0_dim, z0_dim, sigmax_dim, sigmay_dim, sigmaz_dim,
		a0, nwm) = state.namelists.multiwavepackets
	(; rhostrattfc, bvsstrattfc) = state.atmosphere
	#(; u, v) = state.variables.predictands
	#(; zb) = state.wkb

	if random_wavepackets
		construct_random_wavepackets!(state.namelists.multiwavepackets, state.namelists.domain, state.namelists.setting.testcase)
	end

	# Set Coriolis parameter.
	fc = coriolis_frequency * tref

	for iwm in 1:nwm

		amp_wkb = a0[iwm]

		if lambdax_dim[iwm] == 0.0
			wnrk_init = 0.0
		else
			wnrk_init = 2.0 * pi / lambdax_dim[iwm] * lref
		end

		if lambday_dim[iwm] == 0.0
			wnrl_init = 0.0
		else
			wnrl_init = 2.0 * pi / lambday_dim[iwm] * lref
		end

		wnrm_init = 2.0 * pi / lambdaz_dim[iwm] * lref

		wnrh_init = sqrt(wnrk_init ^ 2.0 + wnrl_init ^ 2.0)

		sigwpx = sigmax_dim[iwm] / lref
		sigwpy = sigmay_dim[iwm] / lref
		sigwpz = sigmaz_dim[iwm] / lref

		zr0 = z0_dim[iwm] / lref
		xr0 = x0_dim[iwm] / lref
		yr0 = y0_dim[iwm] / lref

		for kz in k0:k1, jy in j0:j1, ix in i0:i1

			wnk_ini[iwm, ix, jy, kz] = wnrk_init
			wnl_ini[iwm, ix, jy, kz] = wnrl_init
			wnm_ini[iwm, ix, jy, kz] = wnrm_init

			# Compute local stratification.
			n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())

			# intrinsic frequency
			omi_notop = branchr * sqrt((n2r * wnrh_init ^ 2
										+
										fc ^ 2 * wnrm_init ^ 2) /
									   (wnrh_init ^ 2 + wnrm_init ^ 2))

			# wave-action density
			fld_amp = (amp_wkb / wnrm_init) ^ 2 *
					  (wnrh_init ^ 2 + wnrm_init ^ 2) / (2.0 * wnrh_init ^ 2) *
					  omi_notop * rhostrattfc[ix, jy, kz]

			if abs(ztfc[ix, jy, kz] - zr0) < sigwpz
				fld_amp *= 0.5 * (1.0 + cos(pi * (ztfc[ix, jy, kz] - zr0) / sigwpz))

				if sigwpx > 0.0 && abs(x[ix+io] - xr0) < sigwpx
					fld_amp *= 0.5 * (1.0 + cos(pi * (x[ix+io] - xr0) / sigwpx))
				elseif sigwpx > 0.0
					fld_amp = 0.0
				end

				if sigwpy > 0.0 && abs(y[jy+jo] - yr0) < sigwpy
					fld_amp *= 0.5 * (1.0 + cos(pi * (y[jy+jo] - yr0) / sigwpy))
				elseif sigwpy > 0.0
					fld_amp = 0.0
				end

			else
				fld_amp = 0.0
			end

			omi_ini[iwm, ix, jy, kz] = omi_notop
			wad_ini[iwm, ix, jy, kz] = fld_amp

		end
	end # iwm

	return
end

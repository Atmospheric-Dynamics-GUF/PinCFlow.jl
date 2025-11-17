"""
```julia
initialize_rays!(state::State)
```

Complete the initialization of MSGWaM by dispatching to a test-case-specific method.

```julia
initialize_rays!(state::State, test_case::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
initialize_rays!(state::State, test_case::AbstractWKBTestCase)
```

Complete the initialization of MSGWaM for WKB test cases.

In each grid cell, `wave_modes` wave modes are computed, using e.g. `activate_orographic_source!` for mountain waves. For each of these modes, `nrx * nry * nrz * nrk * nrl * nrm` ray volumes are then defined such that they evenly divide the volume one would get for `nrx = nry = nrz = nrk = nrl = nrm = 1` (the parameters are taken from `state.namelists.wkb`). Finally, the maximum group velocities are determined for the corresponding CFL condition that is used in the computation of the time step.

# Arguments

  - `state`: Model state.

  - `test_case`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)
"""
function initialize_rays! end

function initialize_rays!(state::State)
    (; test_case) = state.namelists.setting
    initialize_rays!(state, test_case)
    return
end

function initialize_rays!(state::State, test_case::AbstractTestCase)
    return
end

function initialize_rays!(state::State, test_case::AbstractWKBTestCase)
    (; x_size, y_size, z_size) = state.namelists.domain
    (; test_case) = state.namelists.setting
    (; coriolis_frequency) = state.namelists.atmosphere
    (;
        xrmin,
        xrmax,
        yrmin,
        yrmax,
        nrx,
        nry,
        nrz,
        nrk,
        nrl,
        nrm,
        wave_modes,
        dkr_factor,
        dlr_factor,
        dmr_factor,
        wkb_mode,
        wave_modes,
    ) = state.namelists.wkb
    (; lref, tref) = state.constants
    (; comm, master, nxx, nyy, nzz, io, jo, ko, i0, i1, j0, j1, k0, k1) =
        state.domain
    (; lx, ly, lz, dx, dy, dz, x, y, zc, jac) = state.grid
    (;
        nray_max,
        nray_wrk,
        n_sfc,
        nray,
        rays,
        surface_indices,
        cgx_max,
        cgy_max,
        cgz_max,
    ) = state.wkb
    (; icesetup) = state.namelists.ice

	# Set Coriolis parameter.
	fc = coriolis_frequency * tref

    # Set zonal index bounds.
    if test_case == WKBMountainWave()
        imin = i0
        imax = i1
    else
        imin = max(i0, floor(Int, (xrmin / lref + lx / 2) / dx) + i0 - io)
        imax = min(i1, floor(Int, (xrmax / lref + lx / 2) / dx) + i0 - io)
    end

    # Set meridional index bounds.
    if test_case == WKBMountainWave()
        jmin = j0
        jmax = j1
    else
        jmin = max(j0, floor(Int, (yrmin / lref + ly / 2) / dy) + j0 - jo)
        jmax = min(j1, floor(Int, (yrmax / lref + ly / 2) / dy) + j0 - jo)
    end

    # Set vertical index bounds.
    if test_case == WKBMountainWave() && ko == 0
        kmin = k0 - 1
        kmax = k0 - 1
    else
        kmin = k0
        kmax = k1
    end

    # Initialize local arrays.
    omi_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnk_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnl_ini = zeros(wave_modes, nxx, nyy, nzz)
    wnm_ini = zeros(wave_modes, nxx, nyy, nzz)
    wad_ini = zeros(wave_modes, nxx, nyy, nzz)

    if wkb_mode == SteadyState() && test_case != WKBMountainWave()
        error(
            "Error in initialize_rays!: SteadyState is implemented for WKBMountainWave only!",
        )
    end

    if test_case == WKBMountainWave()
        activate_orographic_source!(
            state,
            omi_ini,
            wnk_ini,
            wnl_ini,
            wnm_ini,
            wad_ini,
        )
    end
	
	if test_case == WKBMultipleWavePackets()
		activate_multiplewavepackets_source!(
			state,
			omi_ini,
			wnk_ini,
			wnl_ini,
			wnm_ini,
			wad_ini,
		)
	end

	#CHANGES 
	# xrmin, yrmin, etc. dimensional

	dk_ini_nd = 0.0
    dl_ini_nd = 0.0
    dm_ini_nd = 0.0

	# Loop over all grid cells with ray volumes.
	@ivy for k in kmin:kmax, j in jmin:jmax, i in imin:imax
		r = 0
		s = 0

        # Loop over all ray volumes within a spatial cell.
        for ix in 1:nrx,
            ik in 1:nrk,
            jy in 1:nry,
            jl in 1:nrl,
            kz in 1:nrz,
            km in 1:nrm,
            alpha in 1:wave_modes

            # Set ray-volume indices.
            if test_case == WKBMountainWave()
                s += 1

				# Set surface indices.
				surface_indices.ixs[s] = ix
				surface_indices.jys[s] = jy
				surface_indices.kzs[s] = kz
				surface_indices.iks[s] = ik
				surface_indices.jls[s] = jl
				surface_indices.kms[s] = km
				surface_indices.alphas[s] = alpha

				# Set surface ray-volume index.
				if wad_ini[alpha, i, j, k] == 0.0
					surface_indices.rs[s, i, j] = -1
					continue
				else
					r += 1
					surface_indices.rs[s, i, j] = r
				end
			else
				r += 1
			end

            # Set ray-volume positions.
            rays.x[r, i, j, k] = (x[io + i] - 0.5 * dx + (ix - 0.5) * dx / nrx)
            rays.y[r, i, j, k] = (y[jo + j] - 0.5 * dy + (jy - 0.5) * dy / nry)
            rays.z[r, i, j, k] = (
                zc[i, j, k] - 0.5 * jac[i, j, k] * dz +
                (kz - 0.5) * jac[i, j, k] * dz / nrz
            )

			xr = rays.x[r, i, j, k]
			yr = rays.y[r, i, j, k]
			zr = rays.z[r, i, j, k]

			# Check if ray volume is too low.
			if zr < -dz
				error(
					"Error in initialize_rays!: Ray volume",
					r,
					"at",
					i,
					j,
					k,
					"is too low!",
				)
			end

            # Compute local stratification.
			n2r = interpolate_stratification(zr, state, N2())

            # Set spatial extents.
            rays.dxray[r, i, j, k] = dx / nrx
            rays.dyray[r, i, j, k] = dy / nry
            rays.dzray[r, i, j, k] = jac[i, j, k] * dz / nrz

			wnk0 = wnk_ini[alpha, i, j, k]
			wnl0 = wnl_ini[alpha, i, j, k]
			wnm0 = wnm_ini[alpha, i, j, k]

            # Ensure correct wavenumber extents.

            if (test_case == WKBMountainWave() || test_case == WKBMultipleWavePackets()) && x_size > 1
                dk_ini_nd = dkr_factor * sqrt(wnk0^2 + wnl0^2)
            end
            if (test_case == WKBMountainWave() || test_case == WKBMultipleWavePackets()) && y_size > 1
                dl_ini_nd = dlr_factor * sqrt(wnk0^2 + wnl0^2)
            end
            if wnm0 == 0.0
                error("Error in WKB: wnm0 = 0!")
            else
                dm_ini_nd = dmr_factor * abs(wnm0)
            end

            # Set ray-volume wavenumbers.
            rays.k[r, i, j, k] =
                (wnk0 - 0.5 * dk_ini_nd + (ik - 0.5) * dk_ini_nd / nrk)
            rays.l[r, i, j, k] =
                (wnl0 - 0.5 * dl_ini_nd + (jl - 0.5) * dl_ini_nd / nrl)
            rays.m[r, i, j, k] =
                (wnm0 - 0.5 * dm_ini_nd + (km - 0.5) * dm_ini_nd / nrm)

            # Set spectral extents.
            rays.dkray[r, i, j, k] = dk_ini_nd / nrk
            rays.dlray[r, i, j, k] = dl_ini_nd / nrl
            rays.dmray[r, i, j, k] = dm_ini_nd / nrm

            # Set spectral volume.
            pspvol = dm_ini_nd
            if x_size > 1
                pspvol = pspvol * dk_ini_nd
            end
            if y_size > 1
                pspvol = pspvol * dl_ini_nd
            end
            
            # Set phase-space wave-action density.
            rays.dens[r, i, j, k] = wad_ini[alpha, i, j, k] / pspvol

			# Interpolate winds to ray-volume position.
			uxr = interpolate_mean_flow(xr, yr, zr, state, U())
			vyr = interpolate_mean_flow(xr, yr, zr, state, V())
			wzr = interpolate_mean_flow(xr, yr, zr, state, W())

			wnrk = rays.k[r, i, j, k]
			wnrl = rays.l[r, i, j, k]
			wnrm = rays.m[r, i, j, k]
			wnrh = sqrt(wnrk^2 + wnrl^2)
			omir = omi_ini[alpha, i, j, k]

			# Compute maximum group velocities.
			cgirx = wnrk * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
			if abs(uxr + cgirx) > abs(cgx_max[])
				cgx_max[] = abs(uxr + cgirx)
			end
			cgiry = wnrl * (n2r - omir^2) / (omir * (wnrh^2 + wnrm^2))
			if abs(vyr + cgiry) > abs(cgy_max[])
				cgy_max[] = abs(vyr + cgiry)
			end
			cgirz = -wnrm * (omir^2 - fc^2) / (omir * (wnrh^2 + wnrm^2))
			if abs(wzr + cgirz) > abs(cgz_max[i, j, k])
				cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(wzr + cgirz))
			end
		end

		# Set ray-volume count.
		nray[i, j, k] = r
		if r > nray_wrk
			error(
				"Error in initialize_rays!: nray",
				[i, j, k],
				" > nray_wrk =",
				nray_wrk,
			)
		end

        # Check if surface ray-volume count is correct.
        if test_case == WKBMountainWave()
            if s != n_sfc
                error(
                    "Error in initialize_rays!: s =",
                    s,
                    "/= n_sfc =",
                    n_sfc,
                    "at (i, j, k) = ",
                    (i, j, k),
                )
            end
        end
    end

	# Compute global ray-volume count.
	@ivy local_sum = sum(nray[imin:imax, jmin:jmax, kmin:kmax])
	global_sum = MPI.Allreduce(local_sum, +, comm)

	# Print information.
	if master
		println("MSGWaM:")
		println("Global ray-volume count: ", global_sum)
		println("Maximum number of ray volumes per cell: ", nray_max)
		println("")
	end

	return
end

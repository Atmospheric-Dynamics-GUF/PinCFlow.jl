function compute_msgwam_ice! end

function compute_msgwam_ice!(state::State)
	icesetup = state.namelists.ice.icesetup
	compute_msgwam_ice!(state, icesetup)
	return
end

function compute_msgwam_ice!(state::State, icesetup::NoIce)
	return
end

function compute_msgwam_ice!(state::State, icesetup::IceOn)
	(; testcase) = state.namelists.setting
	compute_msgwam_ice!(state, testcase)
	return
end

function compute_msgwam_ice!(state::State, testcase::AbstractTestCase)
	return
end

function compute_msgwam_ice!(state::State, testcase::WKBMultipleWavePackets)
	(; domain, grid) = state
	(; sizex, sizey) = state.namelists.domain
	(; coriolis_frequency) = state.namelists.atmosphere
	(; branchr) = state.namelists.wkb
	(; tref, fr2) = state.constants
	(; i0, i1, j0, j1, k0, k1, io, jo) = domain
	(; dx, dy, dz, x, y, ztildetfc, jac) = grid
	(; rhostrattfc, thetastrattfc) = state.atmosphere
	(; nray, rays, integrals) = state.wkb
	(; sgs) = state.ice
	(; nscx, nscy, nscz) = state.namelists.ice
	(; dxsc, dysc, dzsc) = state.ice.subgrid
	(; kappa, ma2) = state.constants

	# Set Coriolis parameter.
	fc = coriolis_frequency * tref

	# set sgs fields to zero
	for field in fieldnames(SgsGW)		
		getfield(sgs, field) .= 0.0
	end

	for kzrv in (k0-1):(k1+1),
		jyrv in (j0-1):(j1+1),
		ixrv in (i0-1):(i1+1)

		for iray in 1:nray[ixrv, jyrv, kzrv]
			if rays.dens[iray, ixrv, jyrv, kzrv] == 0
				continue
			end

			xr = rays.x[iray, ixrv, jyrv, kzrv]
			yr = rays.y[iray, ixrv, jyrv, kzrv]
			zr = rays.z[iray, ixrv, jyrv, kzrv]

			dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
			dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
			dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

			kr = rays.k[iray, ixrv, jyrv, kzrv]
			lr = rays.l[iray, ixrv, jyrv, kzrv]
			mr = rays.m[iray, ixrv, jyrv, kzrv]

			dkr = rays.dkray[iray, ixrv, jyrv, kzrv]
			dlr = rays.dlray[iray, ixrv, jyrv, kzrv]
			dmr = rays.dmray[iray, ixrv, jyrv, kzrv]

			khr = sqrt(kr^2 + lr^2)

			n2r = interpolate_stratification(zr, state, N2())

			omir = 
					branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

			#cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
			#cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
			#cgirz = -mr * (omir^2 - fc^2) / (omir * (khr^2 + mr^2))

			(ixmin, ixmax, jymin, jymax) =
				compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

			for ix in ixmin:ixmax

				for jy in jymin:jymax

					kzmin =
						get_next_half_level(ix, jy, zr - dzr / 2, domain, grid)
					kzmax =
						get_next_half_level(ix, jy, zr + dzr / 2, domain, grid)


					for kz in kzmin:kzmax

						for ii in 1:nscx
							ii2 = (ix - 1) * nscx + ii
							for jj in 1:nscy
								jj2 = (jy - 1) * nscy + jj
								for kk in 1:nscz
									kk2 = (kz - 1) * nscz + kk

									#subcell center
									xsc = x[ix+io] - dx / 2.0 +
										  (ii - 0.5) * dxsc
									ysc = y[jy+jo] - dy / 2.0 +
										  (jj - 0.5) * dysc
									zsc = ztildetfc[ix, jy, kz-1] +
										  (kk - 0.5) * dzsc * jac[ix, jy, kz]

									#interpolate phase/amplitude RV at cell center

									dxx = xsc - xr
									dyy = ysc - yr
									dzz = zsc - zr

									if abs(dxx) <= dxr / 2 &&
									   abs(dyy) <= dyr / 2 &&
									   abs(dzz) <= dzr / 2

										if sizex > 1
											dxi = (
												min(xr + dxr / 2, xsc + dxsc * 0.5) -
												max(xr - dxr / 2, xsc - dxsc * 0.5)
											)
											fcpspx = dkr * dxi / dxsc
										else
											fcpspx = 1.0
										end

										if sizey > 1
											dyi = (
												min((yr + dyr * 0.5), ysc + dysc * 0.5) -
												max((yr - dyr * 0.5), ysc - dysc * 0.5))
											fcpspy = dlr * dyi / dysc
										else
											fcpspy = 1.0
										end

										dzi = (min((zr + dzr * 0.5),
											zsc + dzsc * jac[ix, jy, kz] * 0.5) -
											   max((zr - dzr * 0.5),
											zsc - dzsc * jac[ix, jy, kz] * 0.5))

										fcpspz = dmr * dzi / jac[ix, jy, kz] / dzsc

										fcpswn = fcpspz * fcpspy * fcpspx


										amprw = sqrt(
											abs(omir) * 2.0 * khr^2 /
											(khr^2 + mr^2) /
											rhostrattfc[ix, jy, kz] * fcpswn *
											rays.dens[iray, ixrv, jyrv, kzrv])


										#to be consisten with LES wavepacket simulation
										#there amplitude of b11 is real with sign from vert. wavenumber
										b11 = amprw / abs(omir / n2r) * sign(mr)
										w10 = Complex(0.0, omir / n2r) * b11 # amplitude w

										theta0 = thetastrattfc[ix, jy, kz]
										theta11 = fr2 * theta0 * b11

										pi12 = Complex(0.0, kappa * ma2 *
															(omir * omir - n2r) / n2r / mr /
															thetastrattfc[ix, jy, kz]) * b11

										# phase at cell center
										dphi = (rays.dphi[iray, ixrv, jyrv, kzrv] + kr * dxx + lr * dyy + mr * dzz)


										thetaPrime = real(theta11 * exp(dphi * 1.0im))
										expPrime = real(pi12 * exp(dphi * 1.0im))
										wPrime = real(w10 * exp(dphi * 1.0im))

										#superimpose fields
										sgs.wwp[ii2, jj2, kk2] = sgs.wwp[ii2, jj2, kk2] + wPrime
										sgs.epp[ii2, jj2, kk2] = sgs.epp[ii2, jj2, kk2] + expPrime
										sgs.thp[ii2, jj2, kk2] = sgs.thp[ii2, jj2, kk2] + thetaPrime
									end

								end
							end
						end

					end
				end
			end
		end
	end
end
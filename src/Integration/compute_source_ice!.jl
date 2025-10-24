function compute_source_ice! end

function compute_source_ice!(state::State)
	icesetup = state.namelists.ice.icesetup
	compute_source_ice!(state, icesetup)
	return
end

function compute_source_ice!(state::State, icesetup::AbstractIce)
	@warn "No specific ice source computation implemented for $(typeof(icesetup)). Skipping."
	return
end

function compute_source_ice!(state::State, icesetup::NoIce)
	return
end

function compute_source_ice!(state::State, icesetup::IceOn)
	(; cloudcover) = state.namelists.ice
	compute_source_ice!(state, cloudcover)
	return
end

function compute_source_ice!(state::State, cloudcover::CloudCoverOn)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; nscx, nscy, nscz) = state.namelists.ice
	#(; dxsc, dysc, dzsc) = state.ice.subgrid
	(; n, q, qv) = state.ice.icepredictands
	(; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource, icepredictands, sgstendencies, sgs, sgspredictands, sgsauxiliaries) = state.ice
	#(; iceauxiliaries) = state.ice
	(; kappainv, pref, gamma) = state.constants
	(; press0_dim) = state.namelists.atmosphere

	p0 = press0_dim / pref

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)

		rqv = icepredictands.qv[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhostrattfc[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pstrattfc[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		NIce = icepredictands.n[i, j, k] # N_v = \rho n

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		#*****
		#store large scale values

		pres_ls = pres
		temp_ls = temp

		for ii in 1:nscx
			ii2 = (i - 1) * nscx + ii
			for jj in 1:nscy
				jj2 = (j - 1) * nscy + jj
				for kk in 1:nscz
					kk2 = (k - 1) * nscz + kk

					expPrime = sgs.epp[ii2, jj2, kk2]
					thetaPrime = sgs.thp[ii2, jj2, kk2]
					#compute p', T'
					#pPrime = PStrat(k) * expPrime
					#HOWEVER to be consistent with other implementation here we use
					pPrime = p0 * (exn_p + expPrime) ^ kappainv - pres_ls
					tPrime = thetaPrime * exn_p + expPrime * theta + thetaPrime * expPrime
					#add GW fluctuations to large-scale fields
					pres = pres_ls + pPrime
					temp = temp_ls + tPrime

					psi = psat_ice(temp, iceconstants)

					rqv = sgspredictands.qv[ii2, jj2, kk2]
					NIce = sgspredictands.n[ii2, jj2, kk2]
					sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

					if sice >= iceconstants.S_c

						sice = iceconstants.S_c #set to critical value

						sgstendencies.dn[ii2, jj2, kk2] = dot_n(sice, rhoMean, iceconstants)
					else
						sgstendencies.dn[ii2, jj2, kk2] = 0.0
					end

					dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)
					sgstendencies.dqv[ii2, jj2, kk2] = dqv
					sgstendencies.dq[ii2, jj2, kk2] = -dqv

					#sgsauxiliaries[ii2, jj2, kk2] = sice #full SIce in RT

				end
			end
		end
		#*****


		# dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)

		# icesource.qvsource[i, j, k] = dqv
		# icesource.qsource[i, j, k] = -dqv

		# iceauxiliaries.iaux1[i, j, k] = sice
		# iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
		# iceauxiliaries.iaux3[i, j, k] = dqv
	end

	return
end

function compute_source_ice!(state::State, cloudcover::CloudCoverOn, large_scale_ice::Bool)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; nscx, nscy, nscz) = state.namelists.ice
	#(; dxsc, dysc, dzsc) = state.ice.subgrid
	(; n, q, qv) = state.ice.icepredictands
	(; dz, jac) = state.grid

	(; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource, icepredictands, sgstendencies, sgs, sgspredictands, sgsauxiliaries) = state.ice
	#(; iceauxiliaries) = state.ice
	(; kappainv, pref, Li_hat, gamma) = state.constants
	(; press0_dim) = state.namelists.atmosphere

	p0 = press0_dim / pref

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)

		rqv = icepredictands.qv[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhostrattfc[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pstrattfc[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		NIce = icepredictands.n[i, j, k] # N_v = \rho npxp

		sice_ls = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		#*****
		#store large scale values

		pres_ls = pres
		temp_ls = temp

		if state.namelists.ice.parameterized_nucleation

		dPdz = (kappainv - Li_hat / temp) * (gamma - 1) *
			   (pstrattfc[i, j, k+1] - pstrattfc[i, j, k-1]) / 2.0 / (dz * jac[i, j, k]) / pstrattfc[i, j, k]

		dotPiPrime = w[i, j, k] * dPdz

		end

		for ii in 1:nscx
			ii2 = (i - 1) * nscx + ii
			for jj in 1:nscy
				jj2 = (j - 1) * nscy + jj
				for kk in 1:nscz
					kk2 = (k - 1) * nscz + kk

					expPrime = sgs.epp[ii2, jj2, kk2]
					thetaPrime = sgs.thp[ii2, jj2, kk2]
					#compute p', T'
					#pPrime = PStrat(k) * expPrime
					#HOWEVER to be consistent with other implementation here we use
					pPrime = p0 * (exn_p + expPrime) ^ kappainv - pres_ls
					tPrime = thetaPrime * exn_p + expPrime * theta + thetaPrime * expPrime
					#add GW fluctuations to large-scale fields
					pres = pres_ls + pPrime
					temp = temp_ls + tPrime

					psi = psat_ice(temp, iceconstants)

					rqv = sgspredictands.qv[ii2, jj2, kk2]
					NIce = sgspredictands.n[ii2, jj2, kk2]
					sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

					if sice >= iceconstants.S_c

						sice = iceconstants.S_c #set to critical value

						if state.namelists.ice.parameterized_nucleation

							dotPiPrime = dotPiPrime + sgs.wwp[ii2, jj2, kk2] * dPdz

							FORgw = dotPiPrime * S_c

							n_pre = NIce / rho_full
							n_post = nIce_param_nuc(n_pre, sice, FORgw, temp, pres, psi)

							sgspredictands.n[ii2, jj2, kk2] = n_post * rho_full # N_ice = \rho n
							sgstendencies.dn[ii2, jj2, kk2] = 0.0 # \dot N_ice=0.
						else
							sgstendencies.dn[ii2, jj2, kk2] = dot_n(sice, rhoMean, iceconstants)
						end
					else
						sgstendencies.dn[ii2, jj2, kk2] = 0.0
					end

					dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)
					sgstendencies.dqv[ii2, jj2, kk2] = dqv
					sgstendencies.dq[ii2, jj2, kk2] = -dqv

					#sgsauxiliaries[ii2, jj2, kk2] = sice #full SIce in RT

				end
			end
		end
		#*****


		# dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)

		# icesource.qvsource[i, j, k] = dqv
		# icesource.qsource[i, j, k] = -dqv

		# iceauxiliaries.iaux1[i, j, k] = sice
		# iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
		# iceauxiliaries.iaux3[i, j, k] = dqv
	end

	return
end

function compute_source_ice!(state::State, cloudcover::CloudCoverOff)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; n, q, qv) = state.ice.icepredictands
	(; dz, jac) = state.grid

	(; rhostrattfc, thetastrattfc, bvsstrattfc, pstrattfc) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, pref, gamma, Li_hat) = state.constants
	(; press0_dim) = state.namelists.atmosphere

	p0 = press0_dim / pref

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pstrattfc[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhostrattfc[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pstrattfc[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		NIce = n[i, j, k] # N_v = \rho n

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		if sice >= iceconstants.S_c

			sice = iceconstants.S_c #set to critical value

			if state.namelists.ice.parameterized_nucleation

				dPdz = (kappainv - Li_hat / temp) * (gamma - 1) *
					   (pstrattfc[i, j, k+1] - pstrattfc[i, j, k-1]) / 2.0 / (dz * jac[i, j, k]) / pstrattfc[i, j, k]

				dotPiPrime = w[i, j, k] * dPdz

				#CHANGES***if state.setting.testcase isa WKBMultipleWavePackets

					wPrime = state.ice.sgs.wwp[i, j, k]

					#new GW forcing term: expressed terms of Exner pressure tendency only and
					#include advection of background pressure by GW vertical vel.

					#***dotPiprime = pres/psi * pRef/PsatIceRef * ( kappaInv - Li_hat / temp ) &
					#     * wPrime * gamma_1 * (PStrat(k+1) - PStrat(k-1))/2./dz / PStrat(k) !ONLY if no topography
					#add extra  /exn_p*( dotExpPrime !+ wPrime * ( piStrat(k+1)-piStrat(k-1))/2./dz )

					dotPiPrime = dotPiPrime + wPrime * dPdz

				##**end

				FORgw = dotPiPrime * S_c

				n_pre = NIce / rho_full
				n_post = nIce_param_nuc(n_pre, SIce, FORgw, temp, pres, psi)

				n[i, j, k] = n_post * rho_full # N_ice = \rho n
				icesource.nsource[i, j, k] = 0.0 # \dot N_ice=0.
			else
				icesource.nsource[i, j, k] = dot_n(sice, rhoMean, iceconstants)
			end
		else
			icesource.nsource[i, j, k] = 0.0
		end

		dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)

		icesource.qvsource[i, j, k] = dqv
		icesource.qsource[i, j, k] = -dqv

		iceauxiliaries.iaux1[i, j, k] = sice
		iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
		iceauxiliaries.iaux3[i, j, k] = dqv
	end

	return
end


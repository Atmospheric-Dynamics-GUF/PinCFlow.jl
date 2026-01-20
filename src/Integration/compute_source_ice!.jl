function compute_source_ice! end

function compute_source_ice!(state::State)
	ice_setup = state.namelists.ice.ice_setup
	compute_source_ice!(state, ice_setup)
	return
end

function compute_source_ice!(state::State, ice_setup::AbstractIce)
	@warn "No specific ice source computation implemented for $(typeof(ice_setup)). Skipping."
	return
end

function compute_source_ice!(state::State, ice_setup::NoIce)
	return
end

function compute_source_ice!(state::State, ice_setup::IceOn)
	(; cloudcover) = state.namelists.ice
	compute_source_ice!(state, cloudcover)
	return
end

function compute_source_ice!(state::State, cloudcover::CloudCoverOn)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; nscx, nscy, nscz) = state.namelists.ice
	(; parameterized_nucleation, parameterized_sgs_q) = state.namelists.ice
	#(; dxsc, dysc, dzsc) = state.ice.subgrid
	(; n, q, qv) = state.ice.icepredictands
	(; rhobar, pbar) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource, icepredictands, sgstendencies, sgs, sgspredictands, sgsauxiliaries) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, pref, gamma) = state.constants
	(; ground_pressure) = state.namelists.atmosphere
	(; Li_hat, S_c) = iceconstants
	(; dz, jac) = state.grid

	p0 = ground_pressure / pref

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pbar[i, j, k] / rho_full

		temp = theta * exn_p

		psi_ls = psat_ice(temp, iceconstants)
		pres_ls = pres
		temp_ls = temp

		NIce_ls = n[i, j, k] # N_v = \rho n

		sice_ls = sat_ratio(rqv, pres_ls, psi_ls, rhoMean, iceconstants)

		if parameterized_nucleation

			dPdz = (kappainv - Li_hat / temp) * (gamma - 1) *
				   (pbar[i, j, k+1] - pbar[i, j, k-1]) / 2.0 / (dz * jac[i, j, k]) / pbar[i, j, k]

			dotPiPrime_ls = w[i, j, k] * dPdz

		end

		if parameterized_sgs_q
			nice_mean = 0.0
			nice_cl = 0.0
			numcell_sc = 0.0
			nice_max = n[i, j, k] / rho_full
			clc_pre = iceauxiliaries.clc[i, j, k]
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

					if parameterized_sgs_q
						n_pre = n[i, j, k] / rho_full
						rqv = qv[i, j, k]
					else
						# CHANGES
						#rqv = sgspredictands.qv[ii2, jj2, kk2]
						rqv = qv[i, j, k]
						#NIce = sgspredictands.n[ii2, jj2, kk2]
						NIce = n[i, j, k] / rho_full
						n_pre = NIce / rho_full
					end

					sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

					if sice >= S_c

						sice = S_c #set to critical value

						if parameterized_nucleation

							dotPiPrime = dotPiPrime_ls + sgs.wwp[ii2, jj2, kk2] * dPdz

							FORgw = dotPiPrime * S_c

							if parameterized_sgs_q
								n_pre = nice_max
							end

							n_post = nIce_param_nuc(n_pre, sice, FORgw, temp, pres, psi, iceconstants)

							if parameterized_sgs_q != true
								sgspredictands.n[ii2, jj2, kk2] = n_post * rho_full # N_ice = \rho n
								sgstendencies.dn[ii2, jj2, kk2] = 0.0 # \dot N_ice=0.
							else
								#CHANGES
								n_post = clc_pre * n_post #!!+ (1.0 - clc_pre) * (n_post + n_pre)
								nice_mean += n_post
								#nice_cl += n_post
								numcell_sc += 1

								if n_post > nice_max
									nice_max = n_post
								end
							end

						else
							sgstendencies.dn[ii2, jj2, kk2] = dot_n(sice, rhoMean, iceconstants)
						end
					else
						if parameterized_sgs_q != true
							sgstendencies.dn[ii2, jj2, kk2] = 0.0
						else
							nice_mean += n_pre
						end
					end

					if parameterized_sgs_q != true
						dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)
						sgstendencies.dqv[ii2, jj2, kk2] = dqv
						sgstendencies.dq[ii2, jj2, kk2] = -dqv

						#sgsauxiliaries[ii2, jj2, kk2] = sice #full SIce in RT
					end

				end
			end
		end

		if parameterized_sgs_q

			# store n_max, n_cl, in auxiliaries
			#nice_cl = nice_cl / max(numcell_cl, 1)
			iceauxiliaries.nmax[i, j, k] = nice_max

			clc = numcell_sc / (nscx*nscy*nscz) # cloud cover fraction
			if clc > 0.0
			#	println(i, " ", j, " ", k, " ", clc, " ", numcell_sc, " ", nscx*nscy*nscz)
			end

			#CHANGES 
			iceauxiliaries.clc[i, j, k] = clc_pre + (1.0 - clc_pre) * clc

			iceauxiliaries.ncl[i, j, k] = nice_cl


			# compute large scale n_ice	
			nice_mean = nice_mean / (nscx*nscy*nscz)
			NIce_ls = nice_mean * rho_full
			n[i, j, k] = NIce_ls # N_ice = \rho n
			icesource.nsource[i, j, k] = 0.0 # \dot N_ice=0.
		end

		dqv = dot_qv(sice_ls, NIce_ls, temp_ls, pres_ls, psi_ls, iceconstants)

		icesource.qvsource[i, j, k] = dqv
		icesource.qsource[i, j, k] = -dqv

		iceauxiliaries.iaux1[i, j, k] = sice_ls
		iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
		iceauxiliaries.iaux3[i, j, k] = dqv
	end

	return
end

# function compute_source_ice!(state::State, cloudcover::CloudCoverOn, large_scale_ice::Bool)
# 	(; i0, i1, j0, j1, k0, k1) = state.domain
# 	(; nscx, nscy, nscz) = state.namelists.ice
# 	#(; dxsc, dysc, dzsc) = state.ice.subgrid
# 	(; n, q, qv) = state.ice.icepredictands
# 	(; dz, jac) = state.grid

# 	(; rhobar, thetabar, bvsstrattfc, pbar) = state.atmosphere
# 	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
# 	(; iceconstants) = state.ice
# 	(; icesource, icepredictands, sgstendencies, sgs, sgspredictands, sgsauxiliaries) = state.ice
# 	#(; iceauxiliaries) = state.ice
# 	(; kappainv, pref, gamma) = state.constants
# 	(; Li_hat) = iceconstants
# 	(; ground_pressure) = state.namelists.atmosphere

# 	p0 = ground_pressure / pref

# 	for k in k0:k1, j in j0:j1, i in i0:i1

# 		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
# 		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

# 		rqv = icepredictands.qv[i, j, k]
# 		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
# 		rhoMean = rhobar[i, j, k]
# 		rho_full = rho[i, j, k] + rhoMean
# 		theta = pbar[i, j, k] / rho_full

# 		temp = theta * exn_p

# 		psi = psat_ice(temp, iceconstants)

# 		NIce = icepredictands.n[i, j, k] # N_v = \rho npxp

# 		sice_ls = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

# 		#*****
# 		#store large scale values

# 		pres_ls = pres
# 		temp_ls = temp

# 		if state.namelists.ice.parameterized_nucleation

# 			dPdz = (kappainv - Li_hat / temp) * (gamma - 1) *
# 				   (pbar[i, j, k+1] - pbar[i, j, k-1]) / 2.0 / (dz * jac[i, j, k]) / pbar[i, j, k]

# 			dotPiPrime = w[i, j, k] * dPdz

# 		end

# 		for ii in 1:nscx
# 			ii2 = (i - 1) * nscx + ii
# 			for jj in 1:nscy
# 				jj2 = (j - 1) * nscy + jj
# 				for kk in 1:nscz
# 					kk2 = (k - 1) * nscz + kk

# 					expPrime = sgs.epp[ii2, jj2, kk2]
# 					thetaPrime = sgs.thp[ii2, jj2, kk2]
# 					#compute p', T'
# 					#pPrime = PStrat(k) * expPrime
# 					#HOWEVER to be consistent with other implementation here we use
# 					pPrime = p0 * (exn_p + expPrime) ^ kappainv - pres_ls
# 					tPrime = thetaPrime * exn_p + expPrime * theta + thetaPrime * expPrime
# 					#add GW fluctuations to large-scale fields
# 					pres = pres_ls + pPrime
# 					temp = temp_ls + tPrime

# 					psi = psat_ice(temp, iceconstants)

# 					rqv = sgspredictands.qv[ii2, jj2, kk2]
# 					NIce = sgspredictands.n[ii2, jj2, kk2]
# 					sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

# 					println("to be finished")
# 					exit(1)

# 					dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)
# 					sgstendencies.dqv[ii2, jj2, kk2] = dqv
# 					sgstendencies.dq[ii2, jj2, kk2] = -dqv

# 					#sgsauxiliaries[ii2, jj2, kk2] = sice #full SIce in RT

# 				end
# 			end
# 		end
# 		#*****


# 		# dqv = dot_qv(sice, NIce, temp, pres, psi, iceconstants)

# 		# icesource.qvsource[i, j, k] = dqv
# 		# icesource.qsource[i, j, k] = -dqv

# 		# iceauxiliaries.iaux1[i, j, k] = sice
# 		# iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
# 		# iceauxiliaries.iaux3[i, j, k] = dqv
# 	end

# 	return
# end

function compute_source_ice!(state::State, cloudcover::CloudCoverOff)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; n, q, qv) = state.ice.icepredictands
	(; dz, jac) = state.grid

	(; rhobar, pbar) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, pref, gamma) = state.constants
	(; Li_hat) = iceconstants
	(; ground_pressure) = state.namelists.atmosphere

	p0 = ground_pressure / pref


	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pbar[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		NIce = n[i, j, k] # N_v = \rho n

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		#changes
		iceauxiliaries.iaux1[i, j, k] = sice	

		if sice >= iceconstants.S_c

			sice = iceconstants.S_c #set to critical value

			if state.namelists.ice.parameterized_nucleation

				dPdz = (kappainv - Li_hat / temp) * (gamma - 1) *
					   (pbar[i, j, k+1] - pbar[i, j, k-1]) / 2.0 / (dz * jac[i, j, k]) / pbar[i, j, k]

				dotPiPrime = w[i, j, k] * dPdz

				if state.namelists.ice.ice_test_case isa WKBMultipleWavePackets

					wPrime = state.ice.sgs.wwp[i, j, k]

					#new GW forcing term: expressed terms of Exner pressure tendency only and
					#include advection of background pressure by GW vertical vel.

					#***dotPiprime = pres/psi * pRef/PsatIceRef * ( kappaInv - Li_hat / temp ) &
					#     * wPrime * gamma_1 * (PStrat(k+1) - PStrat(k-1))/2./dz / PStrat(k) !ONLY if no topography
					#add extra  /exn_p*( dotExpPrime !+ wPrime * ( piStrat(k+1)-piStrat(k-1))/2./dz )

					dotPiPrime = dotPiPrime + wPrime * dPdz

				end

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

		
		iceauxiliaries.iaux2[i, j, k] = icesource.nsource[i, j, k]
		iceauxiliaries.iaux3[i, j, k] = dqv
	end

	return
end


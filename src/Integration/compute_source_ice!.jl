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
	(; cloudcover, nucleation) = state.namelists.ice
	compute_source_ice!(state, cloudcover, nucleation)
	return
end

function compute_source_ice!(state::State, cloudcover::AbstractCloudCover, ::AbstractNucleation)
    @warn "No specific nucleation implementation for $(typeof(cloudcover)). Skipping."
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

# combind nucleation: homogeneous and heterogeneous nucleation with separate categories
function compute_source_ice!(state::State, cloudcover::CloudCoverOff, ::BothNucleation)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; n_hom, q_hom, qv, n_in, n_het, q_het) = state.ice.icepredictands
	(; dz, jac, zc) = state.grid

	(; rhobar, pbar, thetabar) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, kappa, pref, gamma, lref) = state.constants
	(; Li_hat, Td_hat, mi0_hat, Tnuc_hat, epsil0, aw, bw_hat, Tr_hat, tRef, p0v_hat, lRef, T0_hat, Ni0_hat, thetaRef, Nref, mRef, rhoRef, mi0, ai, bi_hat) = iceconstants
	(; ground_pressure) = state.namelists.atmosphere
	(; dt_ice) = state.namelists.ice

	p0 = ground_pressure / pref	

	dt_ice = dt_ice / tRef
	
	for k in k0:k1, j in j0:j1, i in i0:i1

		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		rq_het = q_het[i , j, k]
		rq_hom = q_hom[i, j, k]
		N_in = n_in[i, j, k] # N_in = rho * n_in [1/m^3]
		N_het = n_het[i, j, k]
		N_hom = n_hom[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pbar[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		N_hom = n_hom[i, j, k] # N_v = \rho n

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		iceauxiliaries.iaux1[i, j, k] = sice
		
		# ------------------------------------
		#		homogeneous nucleation
		# ------------------------------------
		dn_hom = 0.0
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

				n_pre = N_hom / rho_full
				n_post = nIce_param_nuc(n_pre, sice, FORgw, temp, iceconstants)

				# Maybe dqv is now different because I do not longer directly add n_hom??
				#n_hom[i, j, k] = n_post * rho_full # N_ice = \rho n
				#icesource.n_homsource[i, j, k] = 0.0 # \dot N_ice=0.

				# Convert the direct state update into a tendency
                dn_hom = (n_post * rho_full - N_hom) / dt_ice
			else
				#icesource.n_homsource[i, j, k] = dot_n(sice, rhoMean, iceconstants)
				dn_hom = dot_n(sice, rhoMean, iceconstants)
			end
		else
			icesource.n_homsource[i, j, k] = 0.0
		end

		dqv_hom = dot_qv(sice, N_hom, temp, pres, psi, iceconstants)
		dq_hom = -dqv_hom

		# ------------------------------------
		#      end homogeneous nucleation
		# ------------------------------------

		# ------------------------------------
		#		heterogeneous nucleation
		# ------------------------------------
		# the parameterization for S_nuc is taken from the COSMO documentation (Doms et. al 2019, eq. 5.101)
		S_nuc_cell = 0.0
		ni_actual = 0.0
		psw = psat_water(temp, iceconstants)
		qv_sat_water = qv_sw(pres, psw, iceconstants)
		qv_sat_ice = qv_si(pres, psi, iceconstants)

		# calculate the source for heterogeneous nucleation 
		if (temp < Td_hat && (rqv/rhoMean) >= qv_sat_ice) || ((Td_hat <= temp) && (temp <= Tnuc_hat) && (rqv/rhoMean) >= qv_sat_water)
			# target state: how many IN-particles should be active in total at this temperature
			ni_target = ni(temp, rho_full, iceconstants)
			# reservoir check: how many unused IN particles are left
			ni_actual = max(0.0, min(ni_target, N_in/rhoMean))

			# calculate the source
			S_nuc_cell = rho_full * (ni_actual*mi0_hat) / dt_ice 
		end
		dqv_het = dot_qv(sice, N_het, temp, pres, psi, iceconstants)
		dq_het = -dqv_het
		dn_het = (ni_actual * rhoMean) / dt_ice

		# ------------------------------------
		#		end heterogeneous nucleation
		# ------------------------------------


        
        
		
        # ------------------------------------
        #     Update Sources Safely
        # ------------------------------------

		# update mass sources
        icesource.q_hetsource[i,j,k] = dq_het + S_nuc_cell
        icesource.q_homsource[i, j, k] = dq_hom
        icesource.qvsource[i, j, k]    = -dq_hom - dq_het - S_nuc_cell

        # update the number concentrations
        icesource.n_hetsource[i, j, k] = dn_het
        icesource.n_insource[i, j, k] = -dn_het
		icesource.n_homsource[i, j, k] = dn_hom

		# update the auxiliaries
		iceauxiliaries.iaux2[i, j, k] = icesource.n_homsource[i, j, k]
		iceauxiliaries.iaux3[i, j, k] = icesource.qvsource[i, j, k] 
		iceauxiliaries.iaux4[i, j, k] = icesource.n_insource[i, j, k] 
		iceauxiliaries.iaux5[i, j, k] = icesource.n_hetsource[i, j, k]		
	end
	
end

# only homogeneous nucleation
# homogeneous nucleation and deposition using q_hom, n_hom, and qv
function compute_source_ice!(state::State, cloudcover::CloudCoverOff, ::HomogeneousOnly)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; n_hom, q_hom, qv) = state.ice.icepredictands
	(; dz, jac) = state.grid

	(; rhobar, pbar) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, pref, gamma) = state.constants
	(; Li_hat, Td_hat, mi0_hat, Tnuc_hat, epsil0, aw, bw_hat, Tr_hat, tRef, p0v_hat, lRef, T0_hat, Ni0_hat, thetaRef, Nref, mRef, rhoRef, mi0, ai, bi_hat) = iceconstants
	(; ground_pressure) = state.namelists.atmosphere
	(; dt_ice) = state.namelists.ice

	p0 = ground_pressure / pref

	dt_ice = dt_ice / tRef

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		rq_hom = q_hom[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pbar[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		N_hom = n_hom[i, j, k] # N_v = \rho n

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)

		dn_hom = 0.0

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

				n_pre = N_hom / rho_full
				n_post = nIce_param_nuc(n_pre, sice, FORgw, temp, iceconstants)

				# Maybe dqv is now different because I do not longer directly add n_hom??
				#n_hom[i, j, k] = n_post * rho_full # N_ice = \rho n
				#icesource.n_homsource[i, j, k] = 0.0 # \dot N_ice=0.

				# Convert the direct state update into a tendency
                dn_hom = (n_post * rho_full - N_hom) / dt_ice
			else
				#icesource.n_homsource[i, j, k] = dot_n(sice, rhoMean, iceconstants)
				dn_hom = dot_n(sice, rhoMean, iceconstants)
			end
		else
			icesource.n_homsource[i, j, k] = 0.0
		end
		dqv_hom = dot_qv(sice, N_hom, temp, pres, psi, iceconstants)
		dq_hom = -dqv_hom

		icesource.n_homsource[i, j, k] = dn_hom
		icesource.qvsource[i, j, k] = -dq_hom
		icesource.q_homsource[i, j, k] = dq_hom
		
		iceauxiliaries.iaux2[i, j, k] = icesource.n_homsource[i, j, k]
		iceauxiliaries.iaux3[i, j, k] = dq_hom
	end

	return
end

# only heterogeneous nucleation
# heterogeneous nucleation and deposition using q_het, n_het, n_in and qv
function compute_source_ice!(state::State, cloudcover::CloudCoverOff, ::HeterogeneousOnly)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; qv, n_in, n_het, q_het) = state.ice.icepredictands
	(; dz, jac, zc) = state.grid

	(; rhobar, pbar, thetabar) = state.atmosphere
	(; rho, rhop, u, v, w, pip, p) = state.variables.predictands
	(; iceconstants) = state.ice
	(; icesource) = state.ice
	(; iceauxiliaries) = state.ice
	(; kappainv, kappa, pref, gamma, lref) = state.constants
	(; Li_hat, Td_hat, mi0_hat, Tnuc_hat, epsil0, aw, bw_hat, Tr_hat, tRef, p0v_hat, lRef, T0_hat, Ni0_hat, thetaRef, Nref, mRef, rhoRef, mi0, ai, bi_hat) = iceconstants
	(; ground_pressure) = state.namelists.atmosphere
	(; dt_ice) = state.namelists.ice

	p0 = ground_pressure / pref

	dt_ice = dt_ice / tRef

	for k in k0:k1, j in j0:j1, i in i0:i1

		# Question exn_p = pi(i, j, k) + (pbar[i, j, k] / p0) ^ (gamma - 1)
		exn_p = pip[i, j, k] + (pbar[i, j, k] / p0) ^ (gamma - 1)

		rqv = qv[i, j, k]
		rq_het = q_het[i , j, k]
		N_in = n_in[i, j, k] # N_in = rho * n_in [1/m^3]
		N_het = n_het[i, j, k]
		pres = p0 * exn_p ^ kappainv #kappaInv = c_p/R
		rhoMean = rhobar[i, j, k]
		rho_full = rho[i, j, k] + rhoMean
		theta = pbar[i, j, k] / rho_full

		temp = theta * exn_p

		psi = psat_ice(temp, iceconstants)

		sice = sat_ratio(rqv, pres, psi, rhoMean, iceconstants)


		# ------------------------------------
		#		heterogeneous nucleation
		# ------------------------------------
		# the parameterization for S_nuc is taken from the COSMO documentation (Doms et. al 2019, eq. 5.101)
		S_nuc_cell = 0.0
		ni_actual = 0.0
		psw = psat_water(temp, iceconstants)
		qv_sat_water = qv_sw(pres, psw, iceconstants)
		qv_sat_ice = qv_si(pres, psi, iceconstants)

		# calculate the source for heterogeneous nucleation 
		if (temp < Td_hat && (rqv/rhoMean) >= qv_sat_ice) || ((Td_hat <= temp) && (temp <= Tnuc_hat) && (rqv/rhoMean) >= qv_sat_water)
			# target state: how many IN-particles should be active in total at this temperature
			ni_target = ni(temp, rho_full, iceconstants)
			# reservoir check: how many unused IN particles are left
			ni_actual = max(0.0, min(ni_target, N_in/rhoMean))

			# calculate the source
			S_nuc_cell = rho_full * (ni_actual*mi0_hat) / dt_ice 
		end

		dqv_het = dot_qv(sice, N_het, temp, pres, psi, iceconstants)
		dq_het = -dqv_het
		dn_het = (ni_actual * rhoMean) / dt_ice

		# ------------------------------------
		#		end heterogeneous nucleation
		# ------------------------------------
		
        # ------------------------------------
        #     Update Sources Safely
        # ------------------------------------

		# update mass sources
        icesource.q_hetsource[i,j,k] = dq_het + S_nuc_cell
        icesource.qvsource[i, j, k]    = - dq_het - S_nuc_cell

        # update the heterogeneous number concentrations
        icesource.n_hetsource[i, j, k] = dn_het
        icesource.n_insource[i, j, k] = -dn_het

		iceauxiliaries.iaux3[i, j, k] = sice
		iceauxiliaries.iaux4[i, j, k] = icesource.n_insource[i, j, k] 
		iceauxiliaries.iaux5[i, j, k] = icesource.n_hetsource[i, j, k] 
		iceauxiliaries.iaux6[i, j, k] = icesource.q_hetsource[i, j, k]
	end
end

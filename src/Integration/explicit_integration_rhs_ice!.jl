
function explicit_integration! end

function explicit_integration_rhs_ice!(
	state::State,
	dtstage::AbstractFloat,
)
	icesetup = state.namelists.ice.icesetup
	return explicit_integration_rhs_ice!(state, dtstage, icesetup)
end

function explicit_integration_rhs_ice!(state::State, dtstage::AbstractFloat, icesetup::NoIce)
	return
end

function explicit_integration_rhs_ice!(
	state::State,
	dtstage::AbstractFloat,
	icesetup::IceOn,
)
	(; nstages) = state.time
	(; tracersetup) = state.namelists.tracer
	(; icesetup, dt_ice) = state.namelists.ice
	(; tref) = state.constants
	(; cloudcover) = state.namelists.ice

	(; constant_advection, hor_adv_vel, dt_ice) = state.namelists.ice
	(; nbx, nby) = state.namelists.domain
	(; parameterized_nucleation, constant_advection, hor_adv_vel, parameterized_sgs_q) = state.namelists.ice
	(; icepredictands) = state.ice
	(; sizexx, sizeyy, sizezz) = state.domain
	(; tref, lref) = state.constants

	n_step_ice = ceil(Int, dtstage * tref / dt_ice)
	dtt_ice = dtstage / n_step_ice

	# do advection with constant velocity
	if constant_advection

		hor_adv_vel = hor_adv_vel .* (tref / lref)  # non-dimensionalize velocity 

		di_adv = floor(Int, hor_adv_vel[1] * dtstage / state.grid.dx)
		dj_adv = floor(Int, hor_adv_vel[2] * dtstage / state.grid.dy)

		di2 = floor(Int, hor_adv_vel[1] * dtstage / state.ice.subgrid.dxsc)
		dj2 = floor(Int, hor_adv_vel[2] * dtstage / state.ice.subgrid.dysc)

		# if di_adv > nbx || dj_adv > nby
		# 	@warn "Large constant advection velocity detected: problem ghost cells and boundary conditions. Consider reducing hor_adv_vel or dtstage."
		# 	exit(1)
		# end

		if state.namelists.domain.npx != 1
			@warn "Constant advection of ice predictands not implemented for MPI yet."
			exit(1)
		end

		for (fd, field) in enumerate(fieldnames(IcePredictands))

			if parameterized_nucleation && fd==1 
				continue
			end

			println("Advecting ice predictand field: ", fd, " by (", di_adv, ", ", dj_adv, ") grid cells.", nbx)

			matrix = circshift(getfield(icepredictands, fd)[:, :, :], (di_adv, 0, 0))
			getfield(icepredictands, fd)[:, :, :] .= matrix
			# apply BC after shift if MPI is used
			# set_boundaries!(state, BoundaryPredictands())

			if state.namelists.ice.cloudcover isa CloudCoverOn && parameterized_sgs_q == false
				# also advect subgrid fields
				matrix_sgs = circshift(getfield(state.ice.sgspredictands, fd)[:, :, :], (di2, 0, 0))
				getfield(state.ice.sgspredictands, fd)[:, :, :] .= matrix_sgs
				#CHANGES:
				#NB:  BC of SGS ICE misisng if MPI used
				# apply BC after shift
				#set_boundaries!(state.ice.subgrid, BoundaryPredictands())
			end


		end

	end

	for ii in 1:n_step_ice
		for rkstage in 1:nstages

			#reconstruct!(state)
			#set_boundaries!(state, BoundaryReconstructions())

			#compute_fluxes!(state, p0)

			#set_boundaries!(state, BoundaryFluxes())

			#save_backups!(state, :rho)

			compute_source_ice!(state)
			update!(state, dtt_ice, rkstage, IceUpdatePhy(), cloudcover)

			#apply_unified_sponge!(
			#    state,
			#    stepfrac[rkstage] * dtstage,
			#    time,
			#    tracersetup,
			#)

			set_boundaries!(state, BoundaryPredictands())

		end
	end

	return
end


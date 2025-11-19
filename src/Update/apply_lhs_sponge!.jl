"""
```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::AbstractVariable,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractPredictand,
>>>>>>> cf395edbf2
)
```

Perform an implicit substep to integrate the Rayleigh-damping term that represents the LHS sponge in the prognostic equation for `variable` by dispatching to the appropriate model-specific method.

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::Rho,
	model::Boussinesq,
)
```

Return in Boussinesq mode (constant density).

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::Rho,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Union{PseudoIncompressible, Compressible},
>>>>>>> cf395edbf2
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the continuity equation.

The update is given by

```math
\\rho \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left(\\rho + \\alpha_\\mathrm{R} \\Delta t \\overline{\\rho}\\right),
```

where ``\\alpha_\\mathrm{R}`` is the Rayleigh-damping coefficient computed by [`PinCFlow.Update.compute_sponges!`](@ref) and ``\\Delta t`` is the time step given as input to this method.

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::RhoP,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::Compressible,
>>>>>>> cf395edbf2
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the auxiliary equation in compressible mode.

The update is given by

```math
\\rho' \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left[\\rho' + \\alpha_\\mathrm{R} \\Delta t \\overline{\\rho} \\left(1 - \\frac{P}{\\rho \\overline{\\theta}}\\right)\\right].
```

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::U,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::Union{Boussinesq, PseudoIncompressible},
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the auxiliary equation in non-compressible modes.

The update is given by

```math
\\rho' \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\rho'.
```

```julia
apply_lhs_sponge!(
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::U,
    model::AbstractModel,
>>>>>>> cf395edbf2
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the zonal-momentum equation.

The update is given by

```math
u_{i + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\left(u_{i + 1 / 2} + \\alpha_{\\mathrm{R}, i + 1 / 2} \\Delta t u_{\\mathrm{R}, i + 1 / 2}\\right).
```

If `state.namelists.sponge.relax_to_mean` is `false`, ``u_{\\mathrm{R}, i + 1 / 2}`` is computed with `state.namelist.sponge.relaxed_u`. Otherwise, it is replaced with the average of ``u_{i + 1 / 2}`` across the terrain-following coordinate surface.

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::V,
	model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the meridional-momentum equation.

The update is given by

```math
v_{j + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\left(v_{j + 1 / 2} + \\alpha_{\\mathrm{R}, j + 1 / 2} \\Delta t v_{\\mathrm{R}, j + 1 / 2}\\right).
```

If `state.namelists.sponge.relax_to_mean` is `false`, ``v_{\\mathrm{R}, j + 1 / 2}`` is computed with `state.namelist.sponge.relaxed_v`. Otherwise, it is replaced with the average of ``v_{j + 1 / 2}`` across the terrain-following coordinate surface.

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::W,
	model::AbstractModel,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the transformed-vertical-momentum equation.

The update is given by

```math
\\widehat{w}_{k + 1 / 2} \\rightarrow \\left(1 + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t\\right)^{- 1} \\left(\\widehat{w}_{k + 1 / 2} + \\alpha_{\\mathrm{R}, k + 1 / 2} \\Delta t \\widehat{w}_{\\mathrm{R}, k + 1 / 2}\\right).
```

If `state.namelists.sponge.relax_to_mean` is `false`, ``\\widehat{w}_{\\mathrm{R}, k + 1 / 2}`` is computed with the functions `relaxed_u`, `relaxed_v` and `relaxed_w` in `state.namelists.sponge`. Otherwise, it is replaced with the average of ``\\widehat{w}_{k + 1 / 2}`` across the terrain-following coordinate surface.

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::PiP,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::Union{Boussinesq, PseudoIncompressible},
>>>>>>> cf395edbf2
)
```

Return in non-compressible modes (Exner-pressure fluctuations are only updated in the corrector step).

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::PiP,
	model::Compressible,
)
```

Update the Exner-pressure fluctuations to account for the Rayleigh damping applied to the mass-weighted potential temperature.

The update is given by

```math
\\pi' \\rightarrow \\pi' - \\alpha_\\mathrm{R} \\Delta t P \\frac{\\partial \\pi'}{\\partial P} \\left(1 - \\frac{\\overline{\\rho}}{\\rho}\\right).
```

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::P,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Union{Boussinesq, PseudoIncompressible},
>>>>>>> cf395edbf2
)
```

Return in non-compressible modes (mass-weighted potential temperature is constant in time).

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::P,
	model::Compressible,
)
```

Integrate the Rayleigh-damping term that represents the LHS sponge in the thermodynamic-energy equation.

The update is given by

```math
P \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} P \\left(1 + \\alpha_\\mathrm{R} \\Delta t \\frac{\\overline{\\rho}}{\\rho}\\right).
```

```julia
apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	tracer_setup::NoTracer,
)
```

Return for configurations without tracer transport.

```julia
apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	tracer_setup::AbstractTracer,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::TracerOn,
>>>>>>> cf395edbf2
)
```

Integrate the Rayleigh-damping terms that represent the LHS sponge in the tracer equations.

In each tracer equation, the update is given by

```math
\\left(\\rho \\chi\\right) \\rightarrow \\left(1 + \\alpha_\\mathrm{R} \\Delta t\\right)^{- 1} \\left[\\rho \\chi + \\alpha_\\mathrm{R} \\Delta t \\left(\\rho \\chi\\right)^{\\left(0\\right)}\\right].
```

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `time`: Simulation time.

  - `variable`: Variable to apply Rayleigh damping to.

  - `model`: Dynamic equations.

  - `tracer_setup`: General tracer-transport configuration.
"""
function apply_lhs_sponge! end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::AbstractVariable,
)
	(; model) = state.namelists.setting
	apply_lhs_sponge!(state, dt, time, variable, model)
	return
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::AbstractPredictand,
)
    (; model) = state.namelists.atmosphere
    apply_lhs_sponge!(state, dt, time, variable, model)
    return
>>>>>>> cf395edbf2
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::Rho,
	model::Boussinesq,
)
	return
end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::Rho,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::Rho,
    model::Union{PseudoIncompressible, Compressible},
>>>>>>> cf395edbf2
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; rho) = state.variables.predictands

	rhobg = 0.0
	@ivy for k in k0:k1, j in j0:j1, i in i0:i1
		alpha = alphar[i, j, k]
		rhoold = rho[i, j, k]
		beta = 1.0 / (1.0 + alpha * dt)
		rhonew = (1.0 - beta) * rhobg + beta * rhoold
		rho[i, j, k] = rhonew
	end

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::RhoP,
	model::Compressible,
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; rhobar, thetabar) = state.atmosphere
	(; alphar) = state.sponge
	(; rho, rhop, p) = state.variables.predictands

	@ivy for k in k0:k1, j in j0:j1, i in i0:i1
		rhopbg =
			rhobar[i, j, k] * (
				1.0 -
				p[i, j, k] / thetabar[i, j, k] /
				(rho[i, j, k] + rhobar[i, j, k])
			)
		alpha = alphar[i, j, k]
		rhopold = rhop[i, j, k]
		beta = 1.0 / (1.0 + alpha * dt)
		rhopnew = (1.0 - beta) * rhopbg + beta * rhopold
		rhop[i, j, k] = rhopnew
	end

	return
end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::RhoP,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::RhoP,
    model::Union{Boussinesq, PseudoIncompressible},
>>>>>>> cf395edbf2
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; rhop) = state.variables.predictands

	rhobg = 0.0
	@ivy for k in k0:k1, j in j0:j1, i in i0:i1
		alpha = alphar[i, j, k]
		rhoold = rhop[i, j, k]
		beta = 1.0 / (1.0 + alpha * dt)
		rhonew = (1.0 - beta) * rhobg + beta * rhoold
		rhop[i, j, k] = rhonew
	end

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::U,
	model::AbstractModel,
)
<<<<<<< HEAD
	(; x_size, y_size) = state.namelists.domain
	(;
		relax_to_mean,
		perturbation_period,
		perturbation_amplitude,
		relaxation_wind,
	) = state.namelists.sponge
	(; uref, tref) = state.constants
	(; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar, horizontal_mean) = state.sponge
	(; u) = state.variables.predictands
=======
    (; x_size, y_size) = state.namelists.domain
    (; relax_to_mean, relaxed_u) = state.namelists.sponge
    (; lref, tref, uref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; x, y, zc) = state.grid
    (; alphar, horizontal_mean) = state.sponge
    (; u) = state.variables.predictands
>>>>>>> cf395edbf2

	(ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

<<<<<<< HEAD
	horizontal_mean .= 0.0

	# Determine relaxation wind.
	@ivy if relax_to_mean
		horizontal_mean .=
			sum(a -> a / x_size / y_size, u[ii, jj, kk]; dims = (1, 2))[1, 1, :]
		MPI.Allreduce!(horizontal_mean, +, layer_comm)
	else
		ubg = relaxation_wind[1] / uref
		if perturbation_period > 0.0
			ubg =
				ubg * (
					1.0 +
					perturbation_amplitude *
					sin(2.0 * pi * time / perturbation_period * tref)
				)
		end
	end

	# Update the zonal wind.
	@ivy for k in kk
		if relax_to_mean
			ubg = horizontal_mean[k-k0+1]
		end
		for j in jj, i in ii
			alpha = 0.5 * (alphar[i, j, k] + alphar[i+1, j, k])
			uold = u[i, j, k]
			beta = 1.0 / (1.0 + alpha * dt)
			unew = (1.0 - beta) * ubg + beta * uold
			u[i, j, k] = unew
		end
	end
=======
    # Compute the horizontal mean.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / x_size / y_size, u[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
    end

    # Update the zonal wind.
    @ivy for k in kk, j in jj, i in ii
        xldim = x[i] * lref
        xrdim = x[i + 1] * lref
        ydim = y[j] * lref
        zcdim = zc[i, j, k] * lref
        tdim = time * tref
        dtdim = dt * tref

        ubg =
            relax_to_mean ? horizontal_mean[k - k0 + 1] :
            (
                relaxed_u(xldim, ydim, zcdim, tdim, dtdim) +
                relaxed_u(xrdim, ydim, zcdim, tdim, dtdim)
            ) / uref / 2

        alpha = 0.5 * (alphar[i, j, k] + alphar[i + 1, j, k])
        uold = u[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        unew = (1.0 - beta) * ubg + beta * uold
        u[i, j, k] = unew
    end
>>>>>>> cf395edbf2

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::V,
	model::AbstractModel,
)
<<<<<<< HEAD
	(; x_size, y_size) = state.namelists.domain
	(;
		relax_to_mean,
		perturbation_period,
		perturbation_amplitude,
		relaxation_wind,
	) = state.namelists.sponge
	(; uref, tref) = state.constants
	(; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar, horizontal_mean) = state.sponge
	(; v) = state.variables.predictands
=======
    (; x_size, y_size) = state.namelists.domain
    (; relax_to_mean, relaxed_v) = state.namelists.sponge
    (; lref, tref, uref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; x, y, zc) = state.grid
    (; alphar, horizontal_mean) = state.sponge
    (; v) = state.variables.predictands
>>>>>>> cf395edbf2

	(ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

<<<<<<< HEAD
	horizontal_mean .= 0.0

	# Determine relaxation wind.
	@ivy if relax_to_mean
		horizontal_mean .=
			sum(a -> a / x_size / y_size, v[ii, jj, kk]; dims = (1, 2))[1, 1, :]
		MPI.Allreduce!(horizontal_mean, +, layer_comm)
	else
		vbg = relaxation_wind[2] / uref
		if perturbation_period > 0.0
			vbg =
				vbg * (
					1.0 +
					perturbation_amplitude *
					sin(2.0 * pi * time / perturbation_period * tref)
				)
		end
	end

	# Update the meridional wind.
	@ivy for k in kk
		if relax_to_mean
			vbg = horizontal_mean[k-k0+1]
		end
		for j in jj, i in ii
			alpha = 0.5 * (alphar[i, j, k] + alphar[i, j+1, k])
			vold = v[i, j, k]
			beta = 1.0 / (1.0 + alpha * dt)
			vnew = (1.0 - beta) * vbg + beta * vold
			v[i, j, k] = vnew
		end
	end
=======
    # Compute the horizontal mean.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / x_size / y_size, v[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
    end

    # Update the meridional wind.
    @ivy for k in kk, j in jj, i in ii
        xdim = x[i] * lref
        ybdim = y[j] * lref
        yfdim = y[j + 1] * lref
        zcdim = zc[i, j, k] * lref
        tdim = time * tref
        dtdim = dt * tref

        vbg =
            relax_to_mean ? horizontal_mean[k - k0 + 1] :
            (
                relaxed_v(xdim, ybdim, zcdim, tdim, dtdim) +
                relaxed_v(xdim, yfdim, zcdim, tdim, dtdim)
            ) / uref / 2

        alpha = 0.5 * (alphar[i, j, k] + alphar[i, j + 1, k])
        vold = v[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        vnew = (1.0 - beta) * vbg + beta * vold
        v[i, j, k] = vnew
    end
>>>>>>> cf395edbf2

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::W,
	model::AbstractModel,
)
<<<<<<< HEAD
	(; x_size, y_size) = state.namelists.domain
	(;
		relax_to_mean,
		perturbation_period,
		perturbation_amplitude,
		relaxation_wind,
	) = state.namelists.sponge
	(; uref, tref) = state.constants
	(; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar, horizontal_mean) = state.sponge
	(; w) = state.variables.predictands
	(; jac) = state.grid
=======
    (; x_size, y_size) = state.namelists.domain
    (; relax_to_mean, relaxed_u, relaxed_v, relaxed_w) = state.namelists.sponge
    (; lref, tref, uref) = state.constants
    (; layer_comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; alphar, horizontal_mean) = state.sponge
    (; w) = state.variables.predictands
    (; x, y, zc, jac, met) = state.grid
>>>>>>> cf395edbf2

	(ii, jj, kk) = (i0:i1, j0:j1, k0:k1)

<<<<<<< HEAD
	horizontal_mean .= 0.0

	# Determine relaxation wind.
	@ivy if relax_to_mean
		horizontal_mean .=
			sum(a -> a / x_size / y_size, w[ii, jj, kk]; dims = (1, 2))[1, 1, :]
		MPI.Allreduce!(horizontal_mean, +, layer_comm)
	else
		wbg = relaxation_wind[3] / uref
		if perturbation_period > 0.0
			wbg =
				wbg * (
					1.0 +
					perturbation_amplitude *
					sin(2.0 * pi * time / perturbation_period * tref)
				)
		end
	end

	# Update the vertical wind.
	@ivy for k in kk
		if relax_to_mean
			wbg = horizontal_mean[k-k0+1]
		end
		for j in jj, i in ii
			alpha =
				(
					jac[i, j, k+1] * alphar[i, j, k] +
					jac[i, j, k] * alphar[i, j, k+1]
				) / (jac[i, j, k] + jac[i, j, k+1])
			wold = w[i, j, k]
			beta = 1.0 / (1.0 + alpha * dt)
			wnew = (1.0 - beta) * wbg + beta * wold
			w[i, j, k] = wnew
		end
	end
=======
    # Compute the horizontal mean.
    @ivy if relax_to_mean
        horizontal_mean .=
            sum(a -> a / x_size / y_size, w[ii, jj, kk]; dims = (1, 2))[1, 1, :]
        MPI.Allreduce!(horizontal_mean, +, layer_comm)
    end

    # Update the vertical wind.
    @ivy for k in kk, j in jj, i in ii
        xdim = x[i] * lref
        ydim = y[j] * lref
        zcddim = zc[i, j, k] * lref
        zcudim = zc[i, j, k + 1] * lref
        tdim = time * tref
        dtdim = dt * tref

        wbg =
            relax_to_mean ? horizontal_mean[k - k0 + 1] :
            (
                jac[i, j, k] * (
                    met[i, j, k + 1, 1, 3] *
                    relaxed_u(xdim, ydim, zcudim, tdim, dtdim) +
                    met[i, j, k + 1, 2, 3] *
                    relaxed_v(xdim, ydim, zcudim, tdim, dtdim) +
                    relaxed_w(xdim, ydim, zcudim, tdim, dtdim) /
                    jac[i, j, k + 1]
                ) +
                jac[i, j, k + 1] * (
                    met[i, j, k, 1, 3] *
                    relaxed_u(xdim, ydim, zcddim, tdim, dtdim) +
                    met[i, j, k, 2, 3] *
                    relaxed_v(xdim, ydim, zcddim, tdim, dtdim) +
                    relaxed_w(xdim, ydim, zcddim, tdim, dtdim) / jac[i, j, k]
                )
            ) / (jac[i, j, k] + jac[i, j, k + 1]) / uref

        alpha =
            (
                jac[i, j, k + 1] * alphar[i, j, k] +
                jac[i, j, k] * alphar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        wold = w[i, j, k]
        beta = 1.0 / (1.0 + alpha * dt)
        wnew = (1.0 - beta) * wbg + beta * wold
        w[i, j, k] = wnew
    end
>>>>>>> cf395edbf2

	return
end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::PiP,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::PiP,
    model::Union{Boussinesq, PseudoIncompressible},
>>>>>>> cf395edbf2
)
	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::PiP,
	model::Compressible,
)
	(; gamma, rsp, pref) = state.constants
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; rhobar) = state.atmosphere
	(; rho, pip, p) = state.variables.predictands

	@ivy for k in k0:k1, j in j0:j1, i in i0:i1
		dpdpi =
			1 / (gamma - 1) * (rsp / pref)^(1 - gamma) * p[i, j, k]^(2 - gamma)
		pib =
			rhobar[i, j, k] * p[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]) /
			dpdpi
		alpha = alphar[i, j, k]
		pipold = pip[i, j, k]
		pipnew = pipold - alpha * dt * (p[i, j, k] / dpdpi - pib)
		pip[i, j, k] = pipnew
	end

	return
end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::P,
	model::AbstractModel,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    variable::P,
    model::Union{Boussinesq, PseudoIncompressible},
>>>>>>> cf395edbf2
)
	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	variable::P,
	model::Compressible,
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; rhobar) = state.atmosphere
	(; rho, p) = state.variables.predictands

	@ivy for k in k0:k1, j in j0:j1, i in i0:i1
		pb = rhobar[i, j, k] * p[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
		alpha = alphar[i, j, k]
		pold = p[i, j, k]
		beta = 1 / (1 + alpha * dt)
		pnew = (1 - beta) * pb + beta * pold
		p[i, j, k] = pnew
	end

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	tracer_setup::NoTracer,
)
	return
end

function apply_lhs_sponge!(
<<<<<<< HEAD
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	tracer_setup::AbstractTracer,
=======
    state::State,
    dt::AbstractFloat,
    time::AbstractFloat,
    tracer_setup::TracerOn,
>>>>>>> cf395edbf2
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; tracerpredictands) = state.tracer
	(; initialtracer) = state.tracer.tracerauxiliaries

	@ivy for field in fieldnames(TracerPredictands)
		chi = getfield(tracerpredictands, field)[:, :, :]
		for k in k0:k1, j in j0:j1, i in i0:i1
			alpha = alphar[i, j, k]
			chi_old = chi[i, j, k]
			beta = 1.0 / (1.0 + alpha * dt)
			chi_new = (1.0 - beta) * initialtracer[i, j, k] + beta * chi_old
			chi[i, j, k] = chi_new
		end
	end

	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	icesetup::AbstractIce,
)
	return
end

function apply_lhs_sponge!(
	state::State,
	dt::AbstractFloat,
	time::AbstractFloat,
	icesetup::IceOn,
)
	(; i0, i1, j0, j1, k0, k1) = state.domain
	(; alphar) = state.sponge
	(; icepredictands, iceauxiliaries) = state.ice

	if maximum(alphar) != 0.0
		println("Applying LHS sponge to ice variables not working")
        exit(1)
	end

	# @ivy for field in fieldnames(IcePredictands)
	# 	for k in k0:k1, j in j0:j1, i in i0:i1
	# 		alpha = alphar[i, j, k]
	# 		ice_old = getfield(icepredictands, field)[i, j, k]
	# 		beta = 1.0 / (1.0 + alpha * dt)
	# 		ice_new = (1.0 - beta) * iceauxiliaries[i, j, k] + beta * ice_old
	# 		getfield(icepredictands, field)[i, j, k] = ice_new
	# 	end
	# end

	return
end

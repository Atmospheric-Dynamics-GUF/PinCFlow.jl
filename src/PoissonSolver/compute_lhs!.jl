"""
```julia
compute_lhs!(state::State)::AbstractFloat
```

Compute the scaled left-hand side of the Poisson equation and return a reference tolerance for the convergence criterion by dispatching to a model-specific method.

```julia
compute_lhs!(
    state::State,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
```

Compute the scaled left-hand side of the Poisson equation in pseudo-incompressible/Boussinesq mode and return a reference tolerance for the convergence criterion.

The scaled left-hand side is given by

```math
\\begin{align*}
    b & = \\frac{\\sqrt{\\overline{\\rho}}}{P} \\frac{1}{J c_p} \\left\\{\\frac{1}{\\Delta \\widehat{x}} \\left[\\left(J P\\right)_{i + 1 / 2} u_{i + 1 / 2} - \\left(J P\\right)_{i - 1 / 2} u_{i - 1 / 2}\\right]\\right.\\\\
    & \\qquad \\qquad \\quad + \\frac{1}{\\Delta \\widehat{y}} \\left[\\left(J P\\right)_{j + 1 / 2} v_{j + 1 / 2} - \\left(J P\\right)_{j - 1 / 2} v_{j - 1 / 2}\\right]\\\\
    & \\qquad \\qquad \\quad + \\left.\\frac{1}{\\Delta \\widehat{z}} \\left[\\left(J P\\right)_{k + 1 / 2} \\widehat{w}_{k + 1 / 2} - \\left(J P\\right)_{k - 1 / 2} \\widehat{w}_{k - 1 / 2}\\right]\\right\\}
\\end{align*}
```

and the reference tolerance is given by

```math
\\tau_\\mathrm{ref} = \\frac{\\sum_{i, j, k} b_{i, j, k}^2}{\\sum_{i, j, k} \\left(b_{u, i, j, k}^2 + b_{v, i, j, k}^2 + b_{\\widehat{w}, i, j, k}^2\\right)},
```

where ``b_u``, ``b_v`` and ``b_{\\widehat{w}}`` are the zonal, meridional and vertical parts of ``b``, respectively. Note that in Boussinesq mode, ``P = P_0`` will cancel out, so that the appropriate divergence constraint remains.

```julia
compute_lhs!(state::State, model::Compressible)::AbstractFloat
```

Compute the scaled left-hand side of the Poisson equation in compressible mode and return a reference tolerance for the convergence criterion.

The scaled left-hand side is given by

```math
\\begin{align*}
    b & = \\frac{\\sqrt{\\overline{\\rho}}}{P} \\frac{1}{J c_p} \\left(\\frac{U_{i + 1 / 2} - U_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{V_{j + 1 / 2} - V_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\widehat{W}_{k + 1 / 2} - \\widehat{W}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) - \\frac{\\sqrt{\\overline{\\rho}}}{P} F^P,
\\end{align*}
```

where ``F^P`` is the diabatic heating, as computed by `compute_volume_force`, and the reference tolerance is computed in the same way as in the method for pseudo-incompressible/Boussinesq mode, with ``b_{F, i, j, k}^2`` added to the sum in the denominator, representing the heating term.

# Arguments

  - `state`: Model state.

  - `model`: Dynamic equations.

# See also

  - [`PinCFlow.Update.compute_volume_force`](@ref)
"""
function compute_lhs! end

function compute_lhs!(state::State)::AbstractFloat
    (; model) = state.namelists.setting
    return compute_lhs!(state, model)
end

function compute_lhs!(
    state::State,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
    (; x_size, y_size, z_size) = state.namelists.domain
    (; ma, kappa) = state.constants
    (; comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhobar, pbar) = state.atmosphere
    (; u, v, w) = state.variables.predictands
    (; rhs) = state.poisson

    # Initialize summation variables.
    divl2 = 0.0
    divl2_norm = 0.0

    # Calculate RHS for TFC.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        # Calculate scaling factor.
        fcscal = sqrt(pbar[i, j, k]^2.0 / rhobar[i, j, k])
        # Store velocities at cell edges.
        ur = u[i, j, k]
        ul = u[i - 1, j, k]
        vf = v[i, j, k]
        vb = v[i, j - 1, k]
        wu = w[i, j, k]
        wd = w[i, j, k - 1]
        # Calculate P at cell edges.
        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        pedgel =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i - 1, j, k] * pbar[i - 1, j, k]
            )
        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        pedgeb =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j - 1, k] * pbar[i, j - 1, k]
            )
        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        pedged =
            jac[i, j, k] *
            jac[i, j, k - 1] *
            (pbar[i, j, k] + pbar[i, j, k - 1]) /
            (jac[i, j, k] + jac[i, j, k - 1])
        # Determine indices for RHS.
        ib = i - i0 + 1
        jb = j - j0 + 1
        kb = k - k0 + 1
        # Compute RHS.
        bu = (pedger * ur - pedgel * ul) / dx / jac[i, j, k] * ma^2.0 * kappa
        bv = (pedgef * vf - pedgeb * vb) / dy / jac[i, j, k] * ma^2.0 * kappa
        bw = (pedgeu * wu - pedged * wd) / dz / jac[i, j, k] * ma^2.0 * kappa
        bu /= fcscal
        bv /= fcscal
        bw /= fcscal
        rhs[ib, jb, kb] = bu + bv + bw
        # Compute check sum for solvability criterion.
        divl2 += rhs[ib, jb, kb]^2.0
        divl2_norm += bu^2.0 + bv^2.0 + bw^2.0
    end

    divl2 = MPI.Allreduce(divl2, +, comm)
    divl2 = sqrt(divl2 / x_size / y_size / z_size)

    divl2_norm = MPI.Allreduce(divl2_norm, +, comm)
    divl2_norm = sqrt(divl2_norm / x_size / y_size / z_size)

    if divl2_norm != 0.0
        tolref = divl2 / divl2_norm
    else
        if divl2 == 0.0
            tolref = 1.0
        else
            error("Error in compute_rhs: divl2_norm = 0 while divl2 != 0!")
        end
    end

    return tolref
end

function compute_lhs!(state::State, model::Compressible)::AbstractFloat
    (; x_size, y_size, z_size) = state.namelists.domain
    (; ma, kappa) = state.constants
    (; comm, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhobar, pbar) = state.atmosphere
    (; u, v, w) = state.variables.predictands
    (; rhs) = state.poisson

    # Initialize summation fields.
    divl2 = 0.0
    divl2_norm = 0.0

    # Calculate RHS for TFC.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        # Calculate scaling factor.
        fcscal = sqrt(pbar[i, j, k]^2.0 / rhobar[i, j, k])
        # Store velocities at cell edges.
        ur = u[i, j, k]
        ul = u[i - 1, j, k]
        vf = v[i, j, k]
        vb = v[i, j - 1, k]
        wu = w[i, j, k]
        wd = w[i, j, k - 1]
        # Determine indices for RHS.
        ib = i - i0 + 1
        jb = j - j0 + 1
        kb = k - k0 + 1
        # Compute the heating.
        heating = compute_volume_force(state, i, j, k, P()) * ma^2.0 * kappa
        # Compute RHS.
        bu = (ur - ul) / dx / jac[i, j, k] * ma^2.0 * kappa
        bv = (vf - vb) / dy / jac[i, j, k] * ma^2.0 * kappa
        bw = (wu - wd) / dz / jac[i, j, k] * ma^2.0 * kappa
        bu /= fcscal
        bv /= fcscal
        bw /= fcscal
        heating /= fcscal
        rhs[ib, jb, kb] = bu + bv + bw - heating
        # Compute check sum for solvability criterion.
        divl2 += rhs[ib, jb, kb]^2.0
        divl2_norm += bu^2.0 + bv^2.0 + bw^2.0 + heating^2.0
    end

    divl2 = MPI.Allreduce(divl2, +, comm)
    divl2 = sqrt(divl2 / x_size / y_size / z_size)

    divl2_norm = MPI.Allreduce(divl2_norm, +, comm)
    divl2_norm = sqrt(divl2_norm / x_size / y_size / z_size)

    if divl2_norm != 0.0
        tolref = divl2 / divl2_norm
    else
        if divl2 == 0.0
            tolref = 1.0
        else
            error("Error in compute_rhs: divl2_norm = 0 while divl2 != 0!")
        end
    end

    return tolref
end

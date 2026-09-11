"""
```julia
compute_gw_tracer_integrals!(
    state::State,
    order::Union{LeadingOrder, NextOrder},
    parameters::NamedTuple,
)
```

Compute the leading- or next-order gravity-wave-tracer fluxes by dispatching to the appropriate method.

```julia
compute_gw_tracer_integrals!(
    state::State,
    order::Union{LeadingOrder, NextOrder},
    parameters::NamedTuple,
    tracer_setup::Val{:NoTracer},
)
```

Return for configurations without tracer transport.

```julia
compute_gw_tracer_integrals!(
    state::State,
    order::LeadingOrder,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
```

Compute the leading-order gravity-wave tracer fluxes and wave amplitudes at ``(i, j, k)``.

The zonal, meridional, and vertical fluxes are given by

```math
\\begin{align*}
    \\bar{\\rho}\\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(0)} & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda,\\mu,\\nu} \\left[F \\Re \\left(u^{(0)}_\\mathrm{w}\\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\bar{\\rho}\\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(0)} & = \\frac{\\bar{\\rho}}{2} \\sum_{r,  \\lambda, \\mu, \\nu} \\left[F \\Re \\left(v^{(0)}_\\mathrm{w} \\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\bar{\\rho}\\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(0)} & = \\frac{\\bar{\\rho}}{2} \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(w^{(0)}_\\mathrm{w} \\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},
\\end{align*}
```

with flux contributions given by 

```math 
\\begin{align*}
    \\frac{\\bar{\\rho}}{2} \\Re \\left(u^{(0)}_{\\mathrm{w}, r}\\chi^{{(0)}^*}_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r - m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r\\right], \\\\
    \\frac{\\bar{\\rho}}{2} \\Re \\left(v^{(0)}_{\\mathrm{w}, r} \\chi^{{(0)}^*}_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_{r}} \\frac{m_r}{\\left|\\boldsymbol{k}_{r}\\right|^2} \\mathcal{A}_r \\left[m_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r - k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r\\right], \\\\
    \\frac{\\bar{\\rho}}{2} \\Re \\left(w^{(0)}_{\\mathrm{w}, r} \\chi^{{(0)}^*}_{\\mathrm{w}, r}\\right) & = \\frac{f}{\\hat{\\omega}_r} \\frac{m_r}{\\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r \\left[k_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r - l_r \\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r\\right].
\\end{align*}
```

The leading-order wave amplitudes, used for the next-order flux calculations, are given by 

```math 
\\begin{align*}
    b^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F b^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    u^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F u^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    v^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F v^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    w^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F w^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\pi^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F \\pi^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\chi^{(0)}_\\mathrm{w} & = \\sum_{r, \\lambda,\\mu,\\nu} \\left[F \\chi^{(0)}_{\\mathrm{w}, r}\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},
\\end{align*}
```

with wave contributions

```math 
\\begin{align*}
    b^{(0)}_{\\mathrm{w}, r} & = \\left[\\frac{\\mathcal{A}_r}{\\bar{\\rho}}\\frac{2N_r^2\\left(k_r^2+l_r^2\\right)}{\\hat{\\omega}_r\\left(k_r^2+l_r^2+m_r^2\\right)}\\right]^{1/2}, \\\\
    u^{(0)}_{\\mathrm{w}, r} & = \\frac{i}{m_r N_r^2}\\frac{\\hat{\\omega}_r^2 - N_r^2}{\\hat{\\omega}_r^2 - f^2} \\left(k_r\\hat{\\omega} + i l_r f\\right)b^{(0)}_{\\mathrm{w}, r}, \\\\
    v^{(0)}_{\\mathrm{w}, r} & = \\frac{i}{m_r N_r^2}\\frac{\\hat{\\omega}_r^2 - N_r^2}{\\hat{\\omega}_r^2 - f^2} \\left(l_r\\hat{\\omega} - i k_r f\\right)b^{(0)}_{\\mathrm{w}, r}, \\\\
    w^{(0)}_{\\mathrm{w}, r} & = \\frac{i\\hat{\\omega}}{N_r^2}b^{(0)}_{\\mathrm{w}, r}, \\\\
    \\pi^{(0)}_{\\mathrm{w}, r} & = \\frac{i}{m_r}\\frac{\\hat{\\omega}^2-N_r^2}{N_r^2}\\frac{\\kappa}{\\bar{\\theta}}b^{(0)}_{\\mathrm{w}, r}, \\\\
    \\chi^{(0)}_{\\mathrm{w}, r} & = -\\frac{i}{\\hat{\\omega}}\\left[u^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial x}\\right)_r + v^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial y}\\right)_r + w^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_\\mathrm{b}}{\\partial z}\\right)_r\\right]
\\end{align*}
```

```julia 
function compute_gw_tracer_integrals!(
    state::State,
    order::NextOrder,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
```

Compute the next-order gravity-wave tracer fluxes at ``(i, j, k)``.

The zonal, meridional, and vertical fluxes are given by 

```math
\\begin{align*}
    \\left\\langle \\tilde{u} \\tilde{\\chi} \\right\\rangle^{(1)} & = \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(u^{(0)}_\\mathrm{w}\\chi^{{(1)}^*}_\\mathrm{w} + u^{(1)}_\\mathrm{w}\\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\left\\langle \\tilde{v} \\tilde{\\chi} \\right\\rangle^{(1)} & = \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(v^{(0)}_\\mathrm{w}\\chi^{{(1)}^*}_\\mathrm{w} + v^{(1)}_\\mathrm{w}\\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu},\\\\
    \\left\\langle \\tilde{w} \\tilde{\\chi} \\right\\rangle^{(1)} & = \\sum_{r, \\lambda, \\mu, \\nu} \\left[F \\Re \\left(w^{(0)}_\\mathrm{w}\\chi^{{(1)}^*}_\\mathrm{w} + w^{(1)}_\\mathrm{w}\\chi^{{(0)}^*}_\\mathrm{w}\\right)\\right]_{r, i + \\lambda, j + \\mu, k + \\nu}.
\\end{align*}
```

The next-order wave amplitudes are obtained via the matrix equation 

```math 
\\begin{pmatrix} u^{(1)}_{\\mathrm{w}, r} \\\\ v^{(1)}_{\\mathrm{w}, r} \\\\ w^{(1)}_{\\mathrm{w}, r} \\\\ b^{(1)}_{\\mathrm{w}, r} / N_r \\\\ c_p\\bar{\\theta}\\pi^{(1)}_{\\mathrm{w}, r} \\end{pmatrix} = \\begin{pmatrix} -i\\hat{\\omega}_r & f_0  & 0 & 0 & ik_r \\\\ f_0 & -i\\hat{\\omega}_r & 0 & 0 & il_r \\\\ 0 & 0 & -i\\hat{\\omega} & -N_r & im_r \\\\ 0 & 0 & N_r & -i\\hat{\\omega}_r & 0 \\\\ ik_r & il_r & im_r & 0 & 0 \\end{pmatrix}^{-1}\\begin{pmatrix} R_{u, r} \\\\ R_{v, r} \\\\ R_{w, r} \\\\ R_{b, r} /N_r \\\\ R_{\\pi, r} \\end{pmatrix}
```

The pseudo-inverse of the matrix is calculated using `LinearAlgebra.pinv`. Furthermore, the right-hand sides are given by

```math 
\\begin{align*} 
    R_{u, r} & = -\\left[\\left(\\frac{u^{(0)}_\\mathrm{w} - u^{(0)}_{\\mathrm{w}, \\mathrm{old}}}{\\Delta t}\\right)_r + u_{\\mathrm{b}, r}\\left(\\frac{\\partial u^{(0)}}{\\partial x}\\right)_r + v_{\\mathrm{b}, r}\\left(\\frac{\\partial u^{(0)}}{\\partial y}\\right)_r\\right] \\\\
    & - \\left[u^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right)_r + v^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_r + w^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_r\\right] \\\\
    & - \\kappa \\left[\\bar{\\theta}\\left(\\frac{\\partial \\pi^{(0)}_\\mathrm{w}}{\\partial x}\\right)_r + ik_r \\theta'_r\\pi^{(0)}_{\\mathrm{w}, r} + \\frac{\\bar{\\theta}}{\\hat{g}}\\left(\\frac{\\partial \\pi'}{\\partial x}\\right)_r\\right], \\\\
    R_{v, r} & = -\\left[\\left(\\frac{v^{(0)}_\\mathrm{w} - v^{(0)}_{\\mathrm{w}, \\mathrm{old}}}{\\Delta t}\\right)_r + u_{\\mathrm{b}, r}\\left(\\frac{\\partial v^{(0)}}{\\partial x}\\right)_r + v_{\\mathrm{b}, r}\\left(\\frac{\\partial v^{(0)}}{\\partial y}\\right)_r\\right] \\\\
    & - \\left[u^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_r + v^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right)_r + w^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_r\\right] \\\\
    & - \\kappa \\left[\\bar{\\theta}\\left(\\frac{\\partial \\pi^{(0)}_\\mathrm{w}}{\\partial y}\\right)_r + il_r \\theta'_r\\pi^{(0)}_{\\mathrm{w}, r} + \\frac{\\bar{\\theta}}{\\hat{g}}\\left(\\frac{\\partial \\pi'}{\\partial y}\\right)_r\\right], \\\\
    R_{w, r} & = -\\left[\\left(\\frac{w^{(0)}_\\mathrm{w} - w^{(0)}_{\\mathrm{w}, \\mathrm{old}}}{\\Delta t}\\right)_r + u_{\\mathrm{b}, r}\\left(\\frac{\\partial w^{(0)}}{\\partial x}\\right)_r + v_{\\mathrm{b}, r}\\left(\\frac{\\partial w^{(0)}}{\\partial y}\\right)_r\\right] \\\\
    & - \\kappa \\left[\\bar{\\theta}\\left(\\frac{\\partial \\pi^{(0)}_\\mathrm{w}}{\\partial z}\\right)_r + im_r \\theta'_r\\pi^{(0)}_{\\mathrm{w}, r} + \\frac{\\bar{\\theta}}{\\hat{g}}\\left(\\frac{\\partial \\pi'}{\\partial z}\\right)_r\\right], \\\\
    R_{b, r} & = -\\left[\\left(\\frac{b^{(0)}_\\mathrm{w} - b^{(0)}_{\\mathrm{w}, \\mathrm{old}}}{\\Delta t}\\right)_r + u_{\\mathrm{b}, r}\\left(\\frac{\\partial b^{(0)}}{\\partial x}\\right)_r + v_{\\mathrm{b}, r}\\left(\\frac{\\partial b^{(0)}}{\\partial y}\\right)_r\\right] \\\\
    &- \\frac{1}{\\bar{\\theta}}\\left[u^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\theta'}{\\partial x}\\right)_r + v^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\theta'}{\\partial y}\\right)_r + w^{(0)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\theta'}{\\partial z}\\right)_r\\right], \\\\
    R_{\\pi, r} & = - \\frac{1}{\\bar{\\rho}\\bar{\\theta}}\\left[\\left(\\frac{\\partial \\bar{\\rho}\\bar{\\theta}u^{(0)}_\\mathrm{w}}{\\partial x}\\right)_r + \\left(\\frac{\\partial \\bar{\\rho}\\bar{\\theta}v^{(0)}_\\mathrm{w}}{\\partial y}\\right)_r + \\left(\\frac{\\partial \\bar{\\rho}\\bar{\\theta}w^{(0)}_\\mathrm{w}}{\\partial z}\\right)_r\\right].
\\end{align*}
```

The next-order tracer wave amplitude ``\\chi^{(1)}_{\\mathrm{w}, r}`` is calculated from 

```math 
\\begin{align*}
    \\chi^{(1)}_{\\mathrm{w}, r} = -\\frac{i}{\\hat{\\omega}}\\left[\\left(\\frac{\\chi^{(0)}_\\mathrm{w} - \\chi^{(0)}_{\\mathrm{w}, \\mathrm{old}}}{\\Delta t}\\right)_r + u_{\\mathrm{b}, r}\\left(\\frac{\\partial \\chi^{(0)}}{\\partial x}\\right)_r + v_{\\mathrm{b}, r}\\left(\\frac{\\partial \\chi^{(0)}}{\\partial y}\\right)_r + u^{(1)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_b}{\\partial x}\\right)_r + v^{(1)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_b}{\\partial y}\\right)_r + w^{(1)}_{\\mathrm{w}, r}\\left(\\frac{\\partial \\chi_b}{\\partial z}\\right)_r\\right].
\\end{align*}
```

Here, ``u^{(0)}_{\\mathrm{w}, \\mathrm{old}}``, ``v^{(0)}_{\\mathrm{w}, \\mathrm{old}}``, ``w^{(0)}_{\\mathrm{w}, \\mathrm{old}}``, ``b^{(0)}_{\\mathrm{w}, \\mathrm{old}}``, and ``\\chi^{(0)}_{\\mathrm{w}, \\mathrm{old}}`` are the leading-order wave amplitudes from the previous time-step, stored in `state.tracer.tracerwkbintegrals.uold`, `state.tracer.tracerwkbintegrals.vold`, `state.tracer.tracerwkbintegrals.wold`, `state.tracer.tracerwkbintegrals.bold`, and `state.tracer.tracerwkbintegrals.chiold`, respectively.

# Arguments

  - `state`: Model state.

  - `order`: Computation of the leading- or next-order fluxes.

  - `parameters`: Named tuple containing parameters used for the computations.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_scalar`](@ref)
"""
function compute_gw_tracer_integrals! end

function compute_gw_tracer_integrals!(
    state::State,
    order::Union{LeadingOrder, NextOrder},
    parameters::NamedTuple,
)
    (; tracer_setup) = state.namelists.tracer

    @dispatch_tracer_setup compute_gw_tracer_integrals!(
        state,
        order,
        parameters,
        Val(tracer_setup),
    )
    return
end

function compute_gw_tracer_integrals!(
    state::State,
    order::Union{LeadingOrder, NextOrder},
    parameters::NamedTuple,
    tracer_setup::Val{:NoTracer},
)
    return
end

@ivy function compute_gw_tracer_integrals!(
    state::State,
    order::LeadingOrder,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
    (; uchi0, vchi0, wchi0, uhat, vhat, what, bhat, pihat, chihat) =
        state.tracer.tracerwkbintegrals
    (; leading_order_impact, next_order_impact) = state.namelists.tracer
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; thetabar, rhobar) = state.atmosphere
    (;
        kr,
        lr,
        mr,
        fc,
        omir,
        wadr,
        dens,
        factor,
        dklm,
        xr,
        yr,
        zr,
        iray,
        jray,
        kray,
        n2r,
    ) = parameters
    (; kappa) = state.constants

    if (fc != 0.0 && leading_order_impact) || next_order_impact
        dchidx =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DX())
        dchidy =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DY())
        dchidz =
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DZ())
    end

    if fc != 0.0 && leading_order_impact
        coeff = fc / omir * mr * wadr / (kr^2 + lr^2 + mr^2)

        uchi0[iray, jray, kray] += coeff * (lr * dchidz - mr * dchidy)
        vchi0[iray, jray, kray] += coeff * (mr * dchidx - kr * dchidz)
        wchi0[iray, jray, kray] += coeff * (kr * dchidy - lr * dchidx)
    end

    if next_order_impact
        bamp =
            sqrt(
                dklm * dens / rhobar[iray, jray, kray] *
                2 *
                n2r^2 *
                (kr^2 + lr^2) / (omir * (kr^2 + lr^2 + mr^2)),
            ) / 2
        uamp =
            1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
            (kr * omir + 1im * lr * fc) *
            bamp
        vamp =
            1im / mr / n2r * (omir^2 - n2r) / (omir^2 - fc^2) *
            (lr * omir - 1im * kr * fc) *
            bamp
        wamp = 1im * omir / n2r * bamp
        bhat[iray, jray, kray] += bamp * factor
        uhat[iray, jray, kray] += uamp * factor
        vhat[iray, jray, kray] += vamp * factor
        what[iray, jray, kray] += wamp * factor
        pihat[iray, jray, kray] +=
            1im / mr * (omir^2 - n2r) / n2r / thetabar[iray, jray, kray] *
            kappa *
            bamp *
            factor
        chihat[iray, jray, kray] +=
            -1im / omir *
            (uamp * dchidx + vamp * dchidy + wamp * dchidz) *
            factor
    end

    return
end

@ivy function compute_gw_tracer_integrals!(
    state::State,
    order::NextOrder,
    parameters::NamedTuple,
    tracer_setup::Val{:TracerOn},
)
    (;
        uhat,
        vhat,
        what,
        bhat,
        pihat,
        chihat,
        uhatold,
        vhatold,
        whatold,
        bhatold,
        chihatold,
        uchi1,
        vchi1,
        wchi1,
    ) = state.tracer.tracerwkbintegrals
    (;
        kr,
        lr,
        mr,
        fc,
        omir,
        wadr,
        dens,
        factor,
        xr,
        yr,
        zr,
        iray,
        jray,
        kray,
        n2r,
        dt,
    ) = parameters
    (; rhobar, thetabar, pbar) = state.atmosphere
    (; rho, pip) = state.variables.predictands
    (; chi) = state.tracer.tracerpredictands
    (; kappa, g_ndim) = state.constants

    mat = [
        -1im*omir -fc 0 0 1im*kr
        fc -1im*omir 0 0 1im*lr
        0 0 -1im*omir -sqrt(n2r) 1im*mr
        0 0 sqrt(n2r) -1im*omir 0
        1im*kr 1im*lr 1im*mr 0 0
    ]

    duhatdtr = interpolate_scalar(state, xr, yr, zr, (uhat .- uhatold) ./ dt)
    dvhatdtr = interpolate_scalar(state, xr, yr, zr, (vhat .- vhatold) ./ dt)
    dwhatdtr = interpolate_scalar(state, xr, yr, zr, (what .- whatold) ./ dt)
    dbhatdtr = interpolate_scalar(state, xr, yr, zr, (bhat .- bhatold) ./ dt)
    dchihatdtr =
        interpolate_scalar(state, xr, yr, zr, (chihat .- chihatold) ./ dt)

    uhatr = interpolate_scalar(state, xr, yr, zr, uhat)
    vhatr = interpolate_scalar(state, xr, yr, zr, vhat)
    whatr = interpolate_scalar(state, xr, yr, zr, what)
    bhatr = interpolate_scalar(state, xr, yr, zr, bhat)
    pihatr = interpolate_scalar(state, xr, yr, zr, pihat)
    chihatr = interpolate_scalar(state, xr, yr, zr, chihat)

    ur = interpolate_mean_flow(xr, yr, zr, state, U())
    vr = interpolate_mean_flow(xr, yr, zr, state, V())

    duhatdxr = interpolate_scalar(state, xr, yr, zr, uhat, DX())
    duhatdyr = interpolate_scalar(state, xr, yr, zr, uhat, DY())
    dvhatdxr = interpolate_scalar(state, xr, yr, zr, vhat, DX())
    dvhatdyr = interpolate_scalar(state, xr, yr, zr, vhat, DY())
    dwhatdxr = interpolate_scalar(state, xr, yr, zr, what, DX())
    dwhatdyr = interpolate_scalar(state, xr, yr, zr, what, DY())
    dbhatdxr = interpolate_scalar(state, xr, yr, zr, bhat, DX())
    dbhatdyr = interpolate_scalar(state, xr, yr, zr, bhat, DY())
    dchihatdxr = interpolate_scalar(state, xr, yr, zr, chihat, DX())
    dchihatdyr = interpolate_scalar(state, xr, yr, zr, chihat, DY())

    dudxr = interpolate_mean_flow(xr, yr, zr, state, DUDX())
    dudyr = interpolate_mean_flow(xr, yr, zr, state, DUDY())
    dudzr = interpolate_mean_flow(xr, yr, zr, state, DUDZ())
    dvdxr = interpolate_mean_flow(xr, yr, zr, state, DVDX())
    dvdyr = interpolate_mean_flow(xr, yr, zr, state, DVDY())
    dvdzr = interpolate_mean_flow(xr, yr, zr, state, DVDZ())

    dpihatdxr = interpolate_scalar(state, xr, yr, zr, pihat, DX())
    dpihatdyr = interpolate_scalar(state, xr, yr, zr, pihat, DY())
    dpihatdzr = interpolate_scalar(state, xr, yr, zr, pihat, DZ())

    dpipdxr = interpolate_scalar(state, xr, yr, zr, pip, DX())
    dpipdyr = interpolate_scalar(state, xr, yr, zr, pip, DY())
    dpipdzr = interpolate_scalar(state, xr, yr, zr, pip, DZ())

    theta =
        pbar[iray, jray, kray] /
        (rhobar[iray, jray, kray] + rho[iray, jray, kray]) -
        thetabar[iray, jray, kray]

    dthetadxr = interpolate_scalar(
        state,
        xr,
        yr,
        zr,
        pbar ./ (rhobar .+ rho) .- thetabar,
        DX(),
    )
    dthetadyr = interpolate_scalar(
        state,
        xr,
        yr,
        zr,
        pbar ./ (rhobar .+ rho) .- thetabar,
        DY(),
    )
    dthetadzr = interpolate_scalar(
        state,
        xr,
        yr,
        zr,
        pbar ./ (rhobar .+ rho) .- thetabar,
        DZ(),
    )

    ru =
        -(duhatdtr + ur * duhatdxr + vr * duhatdyr) -
        (uhatr * dudxr + vhatr * dudyr + whatr * dudzr) -
        1 / kappa * (
            thetabar[iray, jray, kray] * dpihatdxr +
            1im * kr * theta * pihatr +
            thetabar[iray, jray, kray] / g_ndim * bhatr * dpipdxr
        )
    rv =
        -(dvhatdtr + ur * dvhatdxr + vr * dvhatdyr) -
        (uhatr * dvdxr + vhatr * dvdyr + whatr * dvdzr) -
        1 / kappa * (
            thetabar[iray, jray, kray] * dpihatdyr +
            1im * lr * theta * pihatr +
            thetabar[iray, jray, kray] / g_ndim * bhatr * dpipdyr
        )
    rw =
        -(dwhatdtr + ur * dwhatdxr + vr * dwhatdyr) -
        1 / kappa * (
            thetabar[iray, jray, kray] * dpihatdzr +
            1im * mr * theta * pihatr +
            thetabar[iray, jray, kray] / g_ndim * bhatr * dpipdzr
        )
    rb =
        -(dbhatdtr + ur * dbhatdxr + vr * dbhatdyr) -
        (uhatr * dthetadxr + vhatr * dthetadyr + whatr * dthetadzr) /
        thetabar[iray, jray, kray]
    rpi =
        -1 / (rhobar[iray, jray, kray] * thetabar[iray, jray, kray]) * (
            interpolate_scalar(
                state,
                xr,
                yr,
                zr,
                rhobar .* thetabar .* uhat,
                DX(),
            ) +
            interpolate_scalar(
                state,
                xr,
                yr,
                zr,
                rhobar .* thetabar .* vhat,
                DY(),
            ) +
            interpolate_scalar(
                state,
                xr,
                yr,
                zr,
                rhobar .* thetabar .* what,
                DZ(),
            )
        )
    matinv = pinv(mat; atol = 1e-12)

    (uhat2r, vhat2r, what2r, bhat2r, pihat2r) =
        matinv * [ru, rv, rw, rb / sqrt(n2r), rpi]

    chihat2r =
        -1im / omir * (
            dchihatdtr +
            ur * dchihatdxr +
            vr * dchihatdyr +
            uhat2r * interpolate_scalar(
                state,
                xr,
                yr,
                zr,
                chi ./ (rho .+ rhobar),
                DX(),
            ) +
            vhat2r * interpolate_scalar(
                state,
                xr,
                yr,
                zr,
                chi ./ (rho .+ rhobar),
                DY(),
            ) +
            what2r *
            interpolate_scalar(state, xr, yr, zr, chi ./ (rho .+ rhobar), DZ())
        )

    uchi1[iray, jray, kray] +=
        real(uhatr * conj(chihat2r) + uhat2r * conj(chihatr)) * factor
    vchi1[iray, jray, kray] +=
        real(vhatr * conj(chihat2r) + vhat2r * conj(chihatr)) * factor
    wchi1[iray, jray, kray] +=
        real(whatr * conj(chihat2r) + what2r * conj(chihatr)) * factor

    return
end

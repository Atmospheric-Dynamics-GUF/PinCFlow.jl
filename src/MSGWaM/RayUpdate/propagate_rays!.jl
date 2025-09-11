"""
```julia
propagate_rays!(state::State, dt::AbstractFloat, rkstage::Integer)
```

Integrate the wave-action-density and ray equations by dispatching to a test-case-specific method.

```julia
propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractTestCase,
)
```

Return for non-WKB test cases.

```julia
propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractWKBTestCase,
)
```

Integrate the wave-action-density and ray equations by dispatching to a WKB-mode-specific method.

```julia
propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::AbstractWKBMode,
)
```

Integrate the wave-action-density and ray equations derived from 1D or 3D transient WKB theory.

The ray equations are given by

```math
\\begin{align*}
    \\frac{\\mathrm{d} x_\\alpha}{\\mathrm{d} t} & = c_{\\mathrm{g}, x, \\alpha} = u_{\\mathrm{b}, \\alpha} + k_\\alpha \\frac{N_\\alpha^2 - \\widehat{\\omega}_\\alpha^2}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\frac{\\mathrm{d} y_\\alpha}{\\mathrm{d} t} & = c_{\\mathrm{g}, y, \\alpha} = v_{\\mathrm{b}, \\alpha} + l_\\alpha \\frac{N_\\alpha^2 - \\widehat{\\omega}_\\alpha^2}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\frac{\\mathrm{d} z_\\alpha}{\\mathrm{d} t} & = c_{\\mathrm{g}, z, \\alpha} = - \\frac{m_\\alpha \\left(\\widehat{\\omega}_\\alpha^2 - f^2\\right)}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2},\\\\
    \\frac{\\mathrm{d} k_\\alpha}{\\mathrm{d} t} & = \\dot{k}_\\alpha = - k_\\alpha \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right)_\\alpha - l_\\alpha \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_\\alpha,\\\\
    \\frac{\\mathrm{d} l_\\alpha}{\\mathrm{d} t} & = \\dot{l}_\\alpha = - k_\\alpha \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_\\alpha - l_\\alpha \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right)_\\alpha,\\\\
    \\frac{\\mathrm{d} m_\\alpha}{\\mathrm{d} t} & = \\dot{m}_\\alpha = - k_\\alpha \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_\\alpha - l_\\alpha \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_\\alpha - \\frac{k_\\alpha^2 + l_\\alpha^2}{2 \\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2} \\left(\\frac{\\partial N^2}{\\partial z}\\right)_\\alpha,
\\end{align*}
```

where the subscript ``\\alpha`` indicates either a ray-volume property or a mean-flow property interpolated to the ray-volume position, via `interpolate_mean_flow` and `interpolate_stratification`. In addition to these, MSGWaM integrates prognostic equations for the ray-volume extents, given by

```math
\\begin{align*}
    \\frac{\\mathrm{d} \\Delta x_\\alpha}{\\mathrm{d} t} & = \\frac{\\mathrm{d} x_{\\alpha, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} x_{\\alpha, -}}{\\mathrm{d} t} = u_{\\mathrm{b}, \\alpha, +} - u_{\\mathrm{b}, \\alpha, -},\\\\
    \\frac{\\mathrm{d} \\Delta y_\\alpha}{\\mathrm{d} t} & = \\frac{\\mathrm{d} y_{\\alpha, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} y_{\\alpha, -}}{\\mathrm{d} t} = v_{\\mathrm{b}, \\alpha, +} - v_{\\mathrm{b}, \\alpha, -},\\\\
    \\frac{\\mathrm{d} \\Delta z_\\alpha}{\\mathrm{d} t} & = \\frac{\\mathrm{d} z_{\\alpha, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} z_{\\alpha, -}}{\\mathrm{d} t} = c_{\\mathrm{g} z, \\alpha, +} - c_{\\mathrm{g} z, \\alpha, -},
\\end{align*}
```

where ``u_{\\mathrm{b}, \\alpha, \\pm}`` is the interpolation of ``u_\\mathrm{b}`` to ``x_{\\alpha, \\pm} = x_\\alpha \\pm \\Delta x_\\alpha / 2`` and ``v_{\\mathrm{b}, \\alpha, \\pm}`` is the equivalent for ``v_\\mathrm{b}`` in ``y``-direction. In the computation of ``c_{\\mathrm{g} z, \\alpha, \\pm}``, the intrinsic frequency and squared buoyancy frequency are interpolated to ``z_{\\alpha, \\pm} = z_\\alpha \\pm \\Delta z_\\alpha / 2``. The update of the spectral ray-volume extents uses the fact that the surfaces in the ``x``-``k``, ``y``-``l`` and ``z``-``m`` subspaces are conserved. Finally, the prognostic equation for the phase-space wave-action density reads

```math
\\frac{\\mathrm{d} \\mathcal{N}_\\alpha}{\\mathrm{d} t} = - 2 \\alpha_{\\mathrm{R}, \\alpha} \\mathcal{N}_\\alpha,
```

where ``\\alpha_{\\mathrm{R}, \\alpha}`` is the interpolation of the Rayleigh-damping coefficient to the ray-volume position, obtained from `interpolate_sponge`. While the ray equations are integrated with the low-storage third-order Runge-Kutta scheme, the phase-space wave-action density is updated with an implicit substep at the end of each Runge-Kutta stage. The group velocities that are calculated for the propagation in physical space are also used to determine the maxima needed for the WKB-CFL condition used in the time-step computation.

```julia
propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::SteadyState,
)
```

Update the vertical wavenumber and wave-action density, using steady-state WKB theory.

In steady-state mode, the ray volumes are stationary in physical space. In mountain-wave simulations, this method first updates the ray volumes in the launch layer by calling `activate_orographic_source!`. Subsequently, it performs a vertical sweep to update all other ray volumes. Therein, the vertical wavenumber is set to

```math
m_\\alpha = - \\sigma \\sqrt{\\frac{\\left(k_\\alpha^2 + l_\\alpha^2\\right) \\left(N_\\alpha^2 - \\widehat{\\omega}_\\alpha^2\\right)}{\\widehat{\\omega}_\\alpha^2 - f^2}},
```

where ``N_\\alpha^2`` is the squared buoyancy frequency interpolated to the ray-volume position (with `interpolate_stratification`) and ``\\widehat{\\omega}_\\alpha = - k_\\alpha u_\\mathrm{b} - l_\\alpha v_\\mathrm{b}`` is the intrinsic frequency (in the case of mountain waves, for which ``\\omega_\\alpha = 0``). The new wave-action-density field is obtained by integrating

```math
\\frac{\\partial \\mathcal{A}_\\alpha}{\\partial z} = - 2 \\alpha_{\\mathrm{R}, \\alpha} \\mathcal{A}_\\alpha - 2 K \\left|\\boldsymbol{k}_\\alpha\\right|^2 \\mathcal{A}_\\alpha,
```

where ``\\alpha_{\\mathrm{R}, \\alpha}`` is the Rayleigh-damping coefficient interpolated to the ray-volume position (using `interpolate_sponge`) and

```math
K = \\left[2 \\sum\\limits_\\alpha \\frac{J \\Delta \\widehat{z}}{c_{\\mathrm{g}, z, \\alpha}} \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right| \\left|\\boldsymbol{k}_\\alpha\\right|\\right)^2 f_\\alpha\\right]^{- 1} \\max \\left[0, \\sum_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right|\\right)^2 f_\\alpha - \\alpha_\\mathrm{S}^2 N^4\\right]
```

is the turbulent viscosity and diffusivity due to wave breaking (see [`PinCFlow.MSGWaM.RayUpdate.apply_saturation_scheme!`](@ref) for more details). After the first right-hand-side term has been integrated with an implicit step, the second term is integrated with the pseudo-time step ``J \\Delta \\widehat{z} / c_{\\mathrm{g} z, \\alpha}``, which corresponds to the substitution ``\\mathcal{A}_\\alpha \\rightarrow \\left(1 - 2 J \\Delta \\widehat{z} / c_{\\mathrm{g} z, \\alpha} K \\left|\\boldsymbol{k}_\\alpha\\right|^2\\right) \\mathcal{A}_\\alpha``. If the domain is parallelized in the vertical, the integration in vertical subdomains is performed sequentially, with one-way communication providing boundary conditions.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `rkstage`: Runge-Kutta-stage index.

  - `testcase`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_position`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_physical_extent`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.get_spectral_extent`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_mean_flow`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_sponge`](@ref)

  - [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.copy_rays!`](@ref)
"""
function propagate_rays! end

function propagate_rays!(state::State, dt::AbstractFloat, rkstage::Integer)
    (; testcase) = state.namelists.setting
    propagate_rays!(state, dt, rkstage, testcase)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractTestCase,
)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    propagate_rays!(state, dt, rkstage, wkb_mode)
    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::AbstractWKBMode,
)
    (; testcase) = state.namelists.setting
    (; branchr, zmin_wkb_dim) = state.namelists.wkb
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; spongelayer) = state.namelists.sponge
    (; lref, tref) = state.constants
    (; nray_max, nray, cgx_max, cgy_max, cgz_max, rays) = state.wkb
    (; dxray, dyray, dzray, dkray, dlray, dmray, ddxray, ddyray, ddzray) =
        state.wkb.increments
    (; alphark, betark, stepfrac, nstages) = state.time
    (; lz, ztildetfc) = state.grid
    (; ko, k0, k1, j0, j1, i0, i1) = state.domain

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    kz0 = testcase == WKBMountainWave() && ko == 0 ? k0 - 1 : k0
    kz1 = k1

    # Initialize WKB increments at the first RK stage.
    @ivy if rkstage == 1
        for kz in kz0:kz1, jy in j0:j1, ix in i0:i1
            for iray in 1:nray[ix, jy, kz]
                dxray[iray, ix, jy, kz] = 0.0
                dyray[iray, ix, jy, kz] = 0.0
                dzray[iray, ix, jy, kz] = 0.0
                dkray[iray, ix, jy, kz] = 0.0
                dlray[iray, ix, jy, kz] = 0.0
                dmray[iray, ix, jy, kz] = 0.0
                ddxray[iray, ix, jy, kz] = 0.0
                ddyray[iray, ix, jy, kz] = 0.0
                ddzray[iray, ix, jy, kz] = 0.0
            end
        end
    end

    cgx_max[] = 0.0
    cgy_max[] = 0.0
    @ivy cgz_max[i0:i1, j0:j1, kz0:kz1] .= 0.0

    @ivy for kz in kz0:kz1, jy in j0:j1, ix in i0:i1
        nskip = 0
        for iray in 1:nray[ix, jy, kz]
            (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            (dxr, dyr, dzr) = get_physical_extent(rays, (iray, ix, jy, kz))
            (axk, ayl, azm) = get_surfaces(rays, (iray, ix, jy, kz))

            xr1 = xr - dxr / 2
            xr2 = xr + dxr / 2
            yr1 = yr - dyr / 2
            yr2 = yr + dyr / 2
            zr1 = zr - dzr / 2
            zr2 = zr + dzr / 2

            khr = sqrt(kr^2 + lr^2)

            # Skip ray volumes that have left the domain.
            if testcase != WKBMountainWave()
                if zr1 < ztildetfc[ix, jy, k0 - 2]
                    nskip += 1
                    continue
                end
            end

            n2r1 = interpolate_stratification(zr1, state, N2())
            n2r = interpolate_stratification(zr, state, N2())
            n2r2 = interpolate_stratification(zr2, state, N2())

            omir1 =
                branchr * sqrt(n2r1 * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            omir =
                branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            omir2 =
                branchr * sqrt(n2r2 * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)

            if any((n2r1, n2r, n2r2) .<= 0)
                error(
                    "Error in propagate_rays!: Interpolated stratification is negative!",
                )
            end

            if khr <= 0
                error(
                    "Error in propagate_rays!: Horizontal wavenumber is negative!",
                )
            end

            # Compute intrinsic zonal group velocity.
            if sizex > 1
                cgirx = kr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            end

            # Compute intrinsic meridional group velocity.
            if sizey > 1
                cgiry = lr * (n2r - omir^2) / (omir * (khr^2 + mr^2))
            end

            # Compute intrinsic vertical group velocities at the vertical edges.
            cgirz1 = -mr * (omir1^2 - fc^2) / (omir1 * (khr^2 + mr^2))
            cgirz2 = -mr * (omir2^2 - fc^2) / (omir2 * (khr^2 + mr^2))

            #-------------------------------
            #      Change of position
            #-------------------------------

            # Update zonal position.

            if sizex > 1 && kz >= k0 && wkb_mode != SingleColumn()
                uxr1 = interpolate_mean_flow(xr1, yr, zr, state, U())
                uxr2 = interpolate_mean_flow(xr2, yr, zr, state, U())

                cgrx1 = cgirx + uxr1
                cgrx2 = cgirx + uxr2

                cgrx = (cgrx1 + cgrx2) / 2

                f = cgrx
                dxray[iray, ix, jy, kz] =
                    dt * f + alphark[rkstage] * dxray[iray, ix, jy, kz]
                rays.x[iray, ix, jy, kz] +=
                    betark[rkstage] * dxray[iray, ix, jy, kz]

                cgx_max[] = max(cgx_max[], abs(cgrx))
            end

            # Update meridional position.

            if sizey > 1 && kz >= k0 && wkb_mode != SingleColumn()
                vyr1 = interpolate_mean_flow(xr, yr1, zr, state, V())
                vyr2 = interpolate_mean_flow(xr, yr2, zr, state, V())

                cgry1 = cgiry + vyr1
                cgry2 = cgiry + vyr2

                cgry = (cgry1 + cgry2) / 2

                f = cgry
                dyray[iray, ix, jy, kz] =
                    dt * f + alphark[rkstage] * dyray[iray, ix, jy, kz]
                rays.y[iray, ix, jy, kz] +=
                    betark[rkstage] * dyray[iray, ix, jy, kz]

                cgy_max[] = max(cgy_max[], abs(cgry))
            end

            # Update vertical position.

            cgrz1 = cgirz1
            cgrz2 = cgirz2

            cgrz = (cgrz1 + cgrz2) / 2

            f = cgrz
            dzray[iray, ix, jy, kz] =
                dt * f + alphark[rkstage] * dzray[iray, ix, jy, kz]
            rays.z[iray, ix, jy, kz] +=
                betark[rkstage] * dzray[iray, ix, jy, kz]

            cgz_max[ix, jy, kz] = max(cgz_max[ix, jy, kz], abs(cgrz))

            # Refraction is only allowed above zmin_wkb_dim / lref.

            if zr > zmin_wkb_dim / lref

                #-------------------------------
                #      Change of wavenumber
                #-------------------------------

                dudxr = interpolate_mean_flow(xr, yr, zr, state, DUDX())
                dudyr = interpolate_mean_flow(xr, yr, zr, state, DUDY())
                dudzr = interpolate_mean_flow(xr, yr, zr, state, DUDZ())

                dvdxr = interpolate_mean_flow(xr, yr, zr, state, DVDX())
                dvdyr = interpolate_mean_flow(xr, yr, zr, state, DVDY())
                dvdzr = interpolate_mean_flow(xr, yr, zr, state, DVDZ())

                dn2dzr = interpolate_stratification(zr, state, DN2DZ())

                dkdt = -dudxr * kr - dvdxr * lr
                dldt = -dudyr * kr - dvdyr * lr
                dmdt =
                    -dudzr * kr - dvdzr * lr -
                    khr^2 * dn2dzr / (2 * omir * (khr^2 + mr^2))

                dkray[iray, ix, jy, kz] =
                    dt * dkdt + alphark[rkstage] * dkray[iray, ix, jy, kz]
                dlray[iray, ix, jy, kz] =
                    dt * dldt + alphark[rkstage] * dlray[iray, ix, jy, kz]
                dmray[iray, ix, jy, kz] =
                    dt * dmdt + alphark[rkstage] * dmray[iray, ix, jy, kz]

                rays.k[iray, ix, jy, kz] +=
                    betark[rkstage] * dkray[iray, ix, jy, kz]
                rays.l[iray, ix, jy, kz] +=
                    betark[rkstage] * dlray[iray, ix, jy, kz]
                rays.m[iray, ix, jy, kz] +=
                    betark[rkstage] * dmray[iray, ix, jy, kz]

                #-------------------------------
                #      Change of extents
                #-------------------------------

                # Update extents in x and k.

                if sizex > 1 && kz >= k0 && wkb_mode != SingleColumn()
                    ddxdt = cgrx2 - cgrx1

                    ddxray[iray, ix, jy, kz] =
                        dt * ddxdt + alphark[rkstage] * ddxray[iray, ix, jy, kz]

                    rays.dxray[iray, ix, jy, kz] +=
                        betark[rkstage] * ddxray[iray, ix, jy, kz]

                    if rays.dxray[iray, ix, jy, kz] <= 0
                        rays.dxray[iray, ix, jy, kz] *= -1
                    end

                    rays.dkray[iray, ix, jy, kz] =
                        axk / rays.dxray[iray, ix, jy, kz]
                end

                # Update extents in y and l.

                if sizey > 1 && kz >= k0 && wkb_mode != SingleColumn()
                    ddydt = cgry2 - cgry1

                    ddyray[iray, ix, jy, kz] =
                        dt * ddydt + alphark[rkstage] * ddyray[iray, ix, jy, kz]

                    rays.dyray[iray, ix, jy, kz] +=
                        betark[rkstage] * ddyray[iray, ix, jy, kz]

                    if rays.dyray[iray, ix, jy, kz] <= 0
                        rays.dyray[iray, ix, jy, kz] *= -1
                    end

                    rays.dlray[iray, ix, jy, kz] =
                        ayl / rays.dyray[iray, ix, jy, kz]
                end

                # Update extents in z and m.

                ddzdt = cgrz2 - cgrz1

                ddzray[iray, ix, jy, kz] =
                    dt * ddzdt + alphark[rkstage] * ddzray[iray, ix, jy, kz]

                rays.dzray[iray, ix, jy, kz] +=
                    betark[rkstage] * ddzray[iray, ix, jy, kz]

                if rays.dzray[iray, ix, jy, kz] <= 0
                    rays.dzray[iray, ix, jy, kz] *= -1
                end

                rays.dmray[iray, ix, jy, kz] =
                    azm / rays.dzray[iray, ix, jy, kz]
            end
        end

        if nskip > 0
            println(
                nskip,
                " out of ",
                nray[ix, jy, kz],
                " ray volumes have been skipped in propagate_rays!!",
            )
            println("")
        end
    end

    #-------------------------------
    #     Change of wave action
    #-------------------------------

    @ivy if spongelayer
        for kz in k0:k1, jy in j0:j1, ix in i0:i1
            for iray in 1:nray[ix, jy, kz]
                (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                betasponge = 1 / (1 + alphasponge * stepfrac[rkstage] * dt)
                rays.dens[iray, ix, jy, kz] *= betasponge
            end
        end
    end

    if testcase == WKBMountainWave()
        activate_orographic_source!(state)
    end

    return
end

function propagate_rays!(
    state::State,
    dt::AbstractFloat,
    rkstage::Integer,
    wkb_mode::SteadyState,
)
    (; sizex, sizey) = state.namelists.domain
    (; testcase) = state.namelists.setting
    (; coriolis_frequency) = state.namelists.atmosphere
    (; spongelayer) = state.namelists.sponge
    (; branchr, lsaturation, alpha_sat) = state.namelists.wkb
    (; stepfrac) = state.time
    (; tref) = state.constants
    (; comm, sizezz, nzz, nx, ny, ko, k0, k1, j0, j1, i0, i1, down, up) =
        state.domain
    (; dx, dy, dz, ztildetfc, ztfc, jac) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; u, v) = state.variables.predictands
    (; nray, rays) = state.wkb

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    if testcase == WKBMountainWave()
        activate_orographic_source!(state)
    end

    @ivy if ko != 0
        nray_down = zeros(Int, nx, ny)
        MPI.Recv!(nray_down, comm; source = down)
        nray[i0:i1, j0:j1, k0 - 1] .= nray_down

        count = maximum(nray[i0:i1, j0:j1, k0 - 1])
        if count > 0
            fields = fieldnames(Rays)
            rays_down = zeros(length(fields), count, nx, ny)
            MPI.Recv!(rays_down, comm; source = down)
            for (index, field) in enumerate(fields)
                getfield(rays, field)[1:count, i0:i1, j0:j1, k0 - 1] .=
                    rays_down[index, :, :, :]
            end
        end
    end

    # Loop over grid cells.
    @ivy for kz in k0:k1, jy in j0:j1, ix in i0:i1

        # Set the ray-volume count.
        nray[ix, jy, kz] = nray[ix, jy, kz - 1]

        # Set up the saturation integrals.
        integral1 = 0.0
        integral2 = 0.0
        m2b2 = 0.0
        m2b2k2 = 0.0

        # Loop over ray volumes.
        for iray in 1:nray[ix, jy, kz]

            # Prepare the ray volume.
            copy_rays!(rays, (iray, ix, jy, kz - 1), (iray, ix, jy, kz))

            # Skip modes with zero wave-action density.
            if rays.dens[iray, ix, jy, kz - 1] == 0
                continue
            end

            # Set the vertical position (and extent).
            rays.z[iray, ix, jy, kz] =
                ztildetfc[ix, jy, kz - 1] +
                (rays.z[iray, ix, jy, kz - 1] - ztildetfc[ix, jy, kz - 2]) /
                jac[ix, jy, kz - 1] * jac[ix, jy, kz]
            rays.dzray[iray, ix, jy, kz] =
                rays.dzray[iray, ix, jy, kz - 1] * jac[ix, jy, kz] /
                jac[ix, jy, kz - 1]

            # Get the horizontal wavenumbers.
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            khr = sqrt(kr^2 + lr^2)

            # Set the reference level.
            kz0 = ko == 0 ? max(k0, kz - 1) : kz - 1

            # Compute the vertical group velocity at the level below.
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz0],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz0] + u[ix - 1, jy, kz0]) / 2 * kr -
                (v[ix, jy, kz0] + v[ix, jy - 1, kz0]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = rays.m[iray, ix, jy, kz0]
                cgirz0 = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end

            # Compute the local vertical wavenumber and vertical group velocity.
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * kr -
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = -branchr * sqrt(khr^2 * (n2r - omir^2) / (omir^2 - fc^2))
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end

            # Set the local vertical wavenumber.
            rays.m[iray, ix, jy, kz] = mr

            # Set the local wave action density.
            if spongelayer
                (xr, yr, zr) = get_physical_position(rays, (iray, ix, jy, kz))
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                rays.dens[iray, ix, jy, kz] =
                    1 / (
                        1 +
                        alphasponge / cgirz *
                        (rays.z[iray, ix, jy, kz] - rays.z[iray, ix, jy, kz0])
                    ) *
                    cgirz0 *
                    rays.dens[iray, ix, jy, kz0] / cgirz
            else
                rays.dens[iray, ix, jy, kz] =
                    cgirz0 * rays.dens[iray, ix, jy, kz0] / cgirz
            end

            # Cycle if the saturation scheme is turned off.
            if !lsaturation
                continue
            end

            # Get the ray volume extents.
            (dxr, dyr, dzr) = get_physical_extent(rays, (iray, ix, jy, kz))
            (dkr, dlr, dmr) = get_spectral_extent(rays, (iray, ix, jy, kz))

            # Compute the phase space factor.
            dzi = min(dzr, jac[ix, jy, kz] * dz)
            facpsp = dzi / jac[ix, jy, kz] / dz * dmr
            if sizex > 1
                dxi = min(dxr, dx)
                facpsp *= dxi / dx * dkr
            end
            if sizey > 1
                dyi = min(dyr, dy)
                facpsp *= dyi / dy * dlr
            end

            # Compute the saturation integrals.
            integral1 = khr^2 * mr^2 / ((khr^2 + mr^2) * omir) * facpsp
            m2b2 +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] *
                integral1 *
                rays.dens[iray, ix, jy, kz]
            integral2 = khr^2 * mr^2 / omir * facpsp
            m2b2k2 +=
                2 * n2r^2 / rhostrattfc[ix, jy, kz] *
                integral2 *
                rays.dens[iray, ix, jy, kz] *
                jac[ix, jy, kz] *
                dz / cgirz
        end

        # Compute the diffusion coefficient
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if m2b2k2 == 0 || m2b2 < alpha_sat^2 * n2r^2
            diffusion = 0.0
        else
            diffusion = (m2b2 - alpha_sat^2 * n2r^2) / (2 * m2b2k2)
        end

        # Reduce the wave action density.
        for iray in 1:nray[ix, jy, kz]
            if !lsaturation
                continue
            end
            if rays.dens[iray, ix, jy, kz] == 0
                continue
            end
            (kr, lr, mr) = get_spectral_position(rays, (iray, ix, jy, kz))
            khr = sqrt(kr^2 + lr^2)
            n2r = interpolate_stratification(
                rays.z[iray, ix, jy, kz],
                state,
                N2(),
            )
            omir =
                -(u[ix, jy, kz] + u[ix - 1, jy, kz]) / 2 * kr -
                (v[ix, jy, kz] + v[ix, jy - 1, kz]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[iray, ix, jy, kz0] = 0.0
                rays.dens[iray, ix, jy, kz] = 0.0
                continue
            end
            rays.dens[iray, ix, jy, kz] *= max(
                0,
                1 -
                jac[ix, jy, kz] * dz / cgirz * 2 * diffusion * (khr^2 + mr^2),
            )
        end
    end

    @ivy if ko + nzz != sizezz
        nray_up = nray[i0:i1, j0:j1, k1]
        MPI.Send(nray_up, comm; dest = up)

        count = maximum(nray[i0:i1, j0:j1, k1])
        if count > 0
            fields = fieldnames(Rays)
            rays_up = zeros(length(fields), count, nx, ny)
        for (index, field) in enumerate(fields)
                rays_up[index, :, :, :] .=
                    getfield(rays, field)[1:count, i0:i1, j0:j1, k1]
            end
            MPI.Send(rays_up, comm; dest = up)
        end
    end

    return
end

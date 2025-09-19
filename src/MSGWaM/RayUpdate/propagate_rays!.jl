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
    \\frac{\\mathrm{d} x_r}{\\mathrm{d} t} & = c_{\\mathrm{g}, x, r} = u_{\\mathrm{b}, r} + k_r \\frac{N_r^2 - \\widehat{\\omega}_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\frac{\\mathrm{d} y_r}{\\mathrm{d} t} & = c_{\\mathrm{g}, y, r} = v_{\\mathrm{b}, r} + l_r \\frac{N_r^2 - \\widehat{\\omega}_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\frac{\\mathrm{d} z_r}{\\mathrm{d} t} & = c_{\\mathrm{g}, z, r} = - \\frac{m_r \\left(\\widehat{\\omega}_r^2 - f^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\frac{\\mathrm{d} k_r}{\\mathrm{d} t} & = \\dot{k}_r = - k_r \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right)_r - l_r \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_r,\\\\
    \\frac{\\mathrm{d} l_r}{\\mathrm{d} t} & = \\dot{l}_r = - k_r \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_r - l_r \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right)_r,\\\\
    \\frac{\\mathrm{d} m_r}{\\mathrm{d} t} & = \\dot{m}_r = - k_r \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_r - l_r \\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_r - \\frac{k_r^2 + l_r^2}{2 \\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\left(\\frac{\\partial N^2}{\\partial z}\\right)_r,
\\end{align*}
```

where the subscript ``r`` indicates either a ray-volume property or a mean-flow property interpolated to the ray-volume position, via `interpolate_mean_flow` and `interpolate_stratification`. In addition to these, MSGWaM integrates prognostic equations for the ray-volume extents, given by

```math
\\begin{align*}
    \\frac{\\mathrm{d} \\Delta x_r}{\\mathrm{d} t} & = \\frac{\\mathrm{d} x_{r, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} x_{r, -}}{\\mathrm{d} t} = u_{\\mathrm{b}, r, +} - u_{\\mathrm{b}, r, -},\\\\
    \\frac{\\mathrm{d} \\Delta y_r}{\\mathrm{d} t} & = \\frac{\\mathrm{d} y_{r, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} y_{r, -}}{\\mathrm{d} t} = v_{\\mathrm{b}, r, +} - v_{\\mathrm{b}, r, -},\\\\
    \\frac{\\mathrm{d} \\Delta z_r}{\\mathrm{d} t} & = \\frac{\\mathrm{d} z_{r, +}}{\\mathrm{d} t} - \\frac{\\mathrm{d} z_{r, -}}{\\mathrm{d} t} = c_{\\mathrm{g} z, r, +} - c_{\\mathrm{g} z, r, -},
\\end{align*}
```

where ``u_{\\mathrm{b}, r, \\pm}`` is the interpolation of ``u_\\mathrm{b}`` to ``x_{r, \\pm} = x_r \\pm \\Delta x_r / 2`` and ``v_{\\mathrm{b}, r, \\pm}`` is the equivalent for ``v_\\mathrm{b}`` in ``y``-direction. In the computation of ``c_{\\mathrm{g} z, r, \\pm}``, the intrinsic frequency and squared buoyancy frequency are interpolated to ``z_{r, \\pm} = z_r \\pm \\Delta z_r / 2``. The update of the spectral ray-volume extents uses the fact that the surfaces in the ``x``-``k``, ``y``-``l`` and ``z``-``m`` subspaces are conserved. Finally, the prognostic equation for the phase-space wave-action density reads

```math
\\frac{\\mathrm{d} \\mathcal{N}_r}{\\mathrm{d} t} = - 2 \\alpha_{\\mathrm{R}, r} \\mathcal{N}_r,
```

where ``\\alpha_{\\mathrm{R}, r}`` is the interpolation of the Rayleigh-damping coefficient to the ray-volume position, obtained from `interpolate_sponge`. While the ray equations are integrated with the low-storage third-order Runge-Kutta scheme, the phase-space wave-action density is updated with an implicit substep at the end of each Runge-Kutta stage. The group velocities that are calculated for the propagation in physical space are also used to determine the maxima needed for the WKB-CFL condition used in the time-step computation.

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
m_r = - \\sigma \\sqrt{\\frac{\\left(k_r^2 + l_r^2\\right) \\left(N_r^2 - \\widehat{\\omega}_r^2\\right)}{\\widehat{\\omega}_r^2 - f^2}},
```

where ``N_r^2`` is the squared buoyancy frequency interpolated to the ray-volume position (with `interpolate_stratification`) and ``\\widehat{\\omega}_r = - k_r u_\\mathrm{b} - l_r v_\\mathrm{b}`` (in the case of mountain waves, for which ``\\omega_r = 0``). The new wave-action-density field is obtained by integrating

```math
\\frac{\\partial}{\\partial z} \\left(c_{\\mathrm{g} z, r} \\mathcal{A}_r\\right) = - 2 \\alpha_{\\mathrm{R}, r} \\mathcal{A}_r - 2 K \\left|\\boldsymbol{k}_r\\right|^2 \\mathcal{A}_r,
```

where ``\\alpha_{\\mathrm{R}, r}`` is the Rayleigh-damping coefficient interpolated to the ray-volume position (using `interpolate_sponge`) and

```math
K = \\left[2 \\sum\\limits_r \\frac{J \\Delta \\widehat{z}}{c_{\\mathrm{g}, z, r}} \\left(m_r \\left|b_{\\mathrm{w}, r}\\right| \\left|\\boldsymbol{k}_r\\right|\\right)^2 f_r\\right]^{- 1} \\max \\left[0, \\sum_r \\left(m_r \\left|b_{\\mathrm{w}, r}\\right|\\right)^2 f_r - \\alpha_\\mathrm{s}^2 N^4\\right]
```

is the turbulent viscosity and diffusivity due to wave breaking (see [`PinCFlow.MSGWaM.RayUpdate.apply_saturation_scheme!`](@ref) for more details). After the first right-hand-side term has been integrated with the implicit step

```math
\\mathcal{A}_r = \\left[1 + \\frac{2 \\alpha_{\\mathrm{R}, r}}{c_{\\mathrm{g} z, r}} \\left(z_r - z_{r, k - 1}\\right)\\right]^{- 1} \\frac{c_{\\mathrm{g} z, r, k - 1}}{c_{\\mathrm{g} z, r}} \\mathcal{A}_{r, k - 1},
```

the second term is integrated with the pseudo-time step ``J \\Delta \\widehat{z} / c_{\\mathrm{g} z, r}``, which corresponds to the substitution ``\\mathcal{A}_r \\rightarrow \\left(1 - 2 J \\Delta \\widehat{z} / c_{\\mathrm{g} z, r} K \\left|\\boldsymbol{k}_r\\right|^2\\right) \\mathcal{A}_r``. If the domain is parallelized in the vertical, the integration in vertical subdomains is performed sequentially, with one-way communication providing boundary conditions.

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

    kmin = testcase == WKBMountainWave() && ko == 0 ? k0 - 1 : k0
    kmax = k1

    # Initialize WKB increments at the first RK stage.
    @ivy if rkstage == 1
        for k in kmin:kmax, j in j0:j1, i in i0:i1
            for r in 1:nray[i, j, k]
                dxray[r, i, j, k] = 0.0
                dyray[r, i, j, k] = 0.0
                dzray[r, i, j, k] = 0.0
                dkray[r, i, j, k] = 0.0
                dlray[r, i, j, k] = 0.0
                dmray[r, i, j, k] = 0.0
                ddxray[r, i, j, k] = 0.0
                ddyray[r, i, j, k] = 0.0
                ddzray[r, i, j, k] = 0.0
            end
        end
    end

    cgx_max[] = 0.0
    cgy_max[] = 0.0
    @ivy cgz_max[i0:i1, j0:j1, kmin:kmax] .= 0.0

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        nskip = 0
        for r in 1:nray[i, j, k]
            (xr, yr, zr) = get_physical_position(rays, r, i, j, k)
            (kr, lr, mr) = get_spectral_position(rays, r, i, j, k)
            (dxr, dyr, dzr) = get_physical_extent(rays, r, i, j, k)
            (axk, ayl, azm) = get_surfaces(rays, r, i, j, k)

            xr1 = xr - dxr / 2
            xr2 = xr + dxr / 2
            yr1 = yr - dyr / 2
            yr2 = yr + dyr / 2
            zr1 = zr - dzr / 2
            zr2 = zr + dzr / 2

            khr = sqrt(kr^2 + lr^2)

            # Skip ray volumes that have left the domain.
            if testcase != WKBMountainWave()
                if zr1 < ztildetfc[i, j, k0 - 2]
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

            if sizex > 1 && k >= k0 && wkb_mode != SingleColumn()
                uxr1 = interpolate_mean_flow(xr1, yr, zr, state, U())
                uxr2 = interpolate_mean_flow(xr2, yr, zr, state, U())

                cgrx1 = cgirx + uxr1
                cgrx2 = cgirx + uxr2

                cgrx = (cgrx1 + cgrx2) / 2

                f = cgrx
                dxray[r, i, j, k] =
                    dt * f + alphark[rkstage] * dxray[r, i, j, k]
                rays.x[r, i, j, k] += betark[rkstage] * dxray[r, i, j, k]

                cgx_max[] = max(cgx_max[], abs(cgrx))
            end

            # Update meridional position.

            if sizey > 1 && k >= k0 && wkb_mode != SingleColumn()
                vyr1 = interpolate_mean_flow(xr, yr1, zr, state, V())
                vyr2 = interpolate_mean_flow(xr, yr2, zr, state, V())

                cgry1 = cgiry + vyr1
                cgry2 = cgiry + vyr2

                cgry = (cgry1 + cgry2) / 2

                f = cgry
                dyray[r, i, j, k] =
                    dt * f + alphark[rkstage] * dyray[r, i, j, k]
                rays.y[r, i, j, k] += betark[rkstage] * dyray[r, i, j, k]

                cgy_max[] = max(cgy_max[], abs(cgry))
            end

            # Update vertical position.

            cgrz1 = cgirz1
            cgrz2 = cgirz2

            cgrz = (cgrz1 + cgrz2) / 2

            f = cgrz
            dzray[r, i, j, k] = dt * f + alphark[rkstage] * dzray[r, i, j, k]
            rays.z[r, i, j, k] += betark[rkstage] * dzray[r, i, j, k]

            cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(cgrz))

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

                dkray[r, i, j, k] =
                    dt * dkdt + alphark[rkstage] * dkray[r, i, j, k]
                dlray[r, i, j, k] =
                    dt * dldt + alphark[rkstage] * dlray[r, i, j, k]
                dmray[r, i, j, k] =
                    dt * dmdt + alphark[rkstage] * dmray[r, i, j, k]

                rays.k[r, i, j, k] += betark[rkstage] * dkray[r, i, j, k]
                rays.l[r, i, j, k] += betark[rkstage] * dlray[r, i, j, k]
                rays.m[r, i, j, k] += betark[rkstage] * dmray[r, i, j, k]

                #-------------------------------
                #      Change of extents
                #-------------------------------

                # Update extents in x and k.

                if sizex > 1 && k >= k0 && wkb_mode != SingleColumn()
                    ddxdt = cgrx2 - cgrx1

                    ddxray[r, i, j, k] =
                        dt * ddxdt + alphark[rkstage] * ddxray[r, i, j, k]

                    rays.dxray[r, i, j, k] +=
                        betark[rkstage] * ddxray[r, i, j, k]

                    if rays.dxray[r, i, j, k] <= 0
                        rays.dxray[r, i, j, k] *= -1
                    end

                    rays.dkray[r, i, j, k] = axk / rays.dxray[r, i, j, k]
                end

                # Update extents in y and l.

                if sizey > 1 && k >= k0 && wkb_mode != SingleColumn()
                    ddydt = cgry2 - cgry1

                    ddyray[r, i, j, k] =
                        dt * ddydt + alphark[rkstage] * ddyray[r, i, j, k]

                    rays.dyray[r, i, j, k] +=
                        betark[rkstage] * ddyray[r, i, j, k]

                    if rays.dyray[r, i, j, k] <= 0
                        rays.dyray[r, i, j, k] *= -1
                    end

                    rays.dlray[r, i, j, k] = ayl / rays.dyray[r, i, j, k]
                end

                # Update extents in z and m.

                ddzdt = cgrz2 - cgrz1

                ddzray[r, i, j, k] =
                    dt * ddzdt + alphark[rkstage] * ddzray[r, i, j, k]

                rays.dzray[r, i, j, k] += betark[rkstage] * ddzray[r, i, j, k]

                if rays.dzray[r, i, j, k] <= 0
                    rays.dzray[r, i, j, k] *= -1
                end

                rays.dmray[r, i, j, k] = azm / rays.dzray[r, i, j, k]
            end
        end

        if nskip > 0
            println(
                nskip,
                " out of ",
                nray[i, j, k],
                " ray volumes have been skipped in propagate_rays!!",
            )
            println("")
        end
    end

    #-------------------------------
    #     Change of wave action
    #-------------------------------

    @ivy if spongelayer
        for k in k0:k1, j in j0:j1, i in i0:i1
            for r in 1:nray[i, j, k]
                (xr, yr, zr) = get_physical_position(rays, r, i, j, k)
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                betasponge = 1 / (1 + alphasponge * stepfrac[rkstage] * dt)
                rays.dens[r, i, j, k] *= betasponge
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

        local_count = maximum(nray[i0:i1, j0:j1, k0 - 1])
        if local_count > 0
            fields = fieldnames(Rays)
            rays_down = zeros(length(fields), local_count, nx, ny)
            MPI.Recv!(rays_down, comm; source = down)
            for (index, field) in enumerate(fields)
                getfield(rays, field)[1:local_count, i0:i1, j0:j1, k0 - 1] .=
                    rays_down[index, :, :, :]
            end
        end
    end

    # Loop over grid cells.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1

        # Set the ray-volume count.
        nray[i, j, k] = nray[i, j, k - 1]

        # Set up the saturation integrals.
        integral1 = 0.0
        integral2 = 0.0
        m2b2 = 0.0
        m2b2k2 = 0.0

        # Loop over ray volumes.
        for r in 1:nray[i, j, k]

            # Prepare the ray volume.
            copy_rays!(rays, r => r, i => i, j => j, k - 1 => k)

            # Skip modes with zero wave-action density.
            if rays.dens[r, i, j, k - 1] == 0
                continue
            end

            # Set the vertical position (and extent).
            rays.z[r, i, j, k] =
                ztildetfc[i, j, k - 1] +
                (rays.z[r, i, j, k - 1] - ztildetfc[i, j, k - 2]) /
                jac[i, j, k - 1] * jac[i, j, k]
            rays.dzray[r, i, j, k] =
                rays.dzray[r, i, j, k - 1] * jac[i, j, k] / jac[i, j, k - 1]

            # Get the horizontal wavenumbers.
            (kr, lr, mr) = get_spectral_position(rays, r, i, j, k)
            khr = sqrt(kr^2 + lr^2)

            # Set the reference level.
            kref = ko == 0 ? max(k0, k - 1) : k - 1

            # Compute the vertical group velocity at the level below.
            n2r = interpolate_stratification(rays.z[r, i, j, kref], state, N2())
            omir =
                -(u[i, j, kref] + u[i - 1, j, kref]) / 2 * kr -
                (v[i, j, kref] + v[i, j - 1, kref]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = rays.m[r, i, j, kref]
                cgirz0 = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[r, i, j, kref] = 0.0
                rays.dens[r, i, j, k] = 0.0
                continue
            end

            # Compute the local vertical wavenumber and vertical group velocity.
            n2r = interpolate_stratification(rays.z[r, i, j, k], state, N2())
            omir =
                -(u[i, j, k] + u[i - 1, j, k]) / 2 * kr -
                (v[i, j, k] + v[i, j - 1, k]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                mr = -branchr * sqrt(khr^2 * (n2r - omir^2) / (omir^2 - fc^2))
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[r, i, j, kref] = 0.0
                rays.dens[r, i, j, k] = 0.0
                continue
            end

            # Set the local vertical wavenumber.
            rays.m[r, i, j, k] = mr

            # Set the local wave action density.
            if spongelayer
                (xr, yr, zr) = get_physical_position(rays, r, i, j, k)
                alphasponge = 2 * interpolate_sponge(xr, yr, zr, state)
                rays.dens[r, i, j, k] =
                    1 / (
                        1 +
                        alphasponge / cgirz *
                        (rays.z[r, i, j, k] - rays.z[r, i, j, kref])
                    ) *
                    cgirz0 *
                    rays.dens[r, i, j, kref] / cgirz
            else
                rays.dens[r, i, j, k] =
                    cgirz0 * rays.dens[r, i, j, kref] / cgirz
            end

            # Cycle if the saturation scheme is turned off.
            if !lsaturation
                continue
            end

            # Get the ray volume extents.
            (dxr, dyr, dzr) = get_physical_extent(rays, r, i, j, k)
            (dkr, dlr, dmr) = get_spectral_extent(rays, r, i, j, k)

            # Compute the phase space factor.
            dzi = min(dzr, jac[i, j, k] * dz)
            facpsp = dzi / jac[i, j, k] / dz * dmr
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
                2 * n2r^2 / rhostrattfc[i, j, k] *
                integral1 *
                rays.dens[r, i, j, k]
            integral2 = khr^2 * mr^2 / omir * facpsp
            m2b2k2 +=
                2 * n2r^2 / rhostrattfc[i, j, k] *
                integral2 *
                rays.dens[r, i, j, k] *
                jac[i, j, k] *
                dz / cgirz
        end

        # Compute the diffusion coefficient
        n2r = interpolate_stratification(ztfc[i, j, k], state, N2())
        if m2b2k2 == 0 || m2b2 < alpha_sat^2 * n2r^2
            diffusion = 0.0
        else
            diffusion = (m2b2 - alpha_sat^2 * n2r^2) / (2 * m2b2k2)
        end

        # Reduce the wave action density.
        for r in 1:nray[i, j, k]
            if !lsaturation
                continue
            end
            if rays.dens[r, i, j, k] == 0
                continue
            end
            (kr, lr, mr) = get_spectral_position(rays, r, i, j, k)
            khr = sqrt(kr^2 + lr^2)
            n2r = interpolate_stratification(rays.z[r, i, j, k], state, N2())
            omir =
                -(u[i, j, k] + u[i - 1, j, k]) / 2 * kr -
                (v[i, j, k] + v[i, j - 1, k]) / 2 * lr
            if fc < branchr * omir < sqrt(n2r)
                cgirz = mr * (fc^2 - n2r) * khr^2 / omir / (khr^2 + mr^2)^2
            else
                rays.dens[r, i, j, kref] = 0.0
                rays.dens[r, i, j, k] = 0.0
                continue
            end
            rays.dens[r, i, j, k] *= max(
                0,
                1 - jac[i, j, k] * dz / cgirz * 2 * diffusion * (khr^2 + mr^2),
            )
        end
    end

    @ivy if ko + nzz != sizezz
        nray_up = nray[i0:i1, j0:j1, k1]
        MPI.Send(nray_up, comm; dest = up)

        local_count = maximum(nray[i0:i1, j0:j1, k1])
        if local_count > 0
            fields = fieldnames(Rays)
            rays_up = zeros(length(fields), local_count, nx, ny)
            for (index, field) in enumerate(fields)
                rays_up[index, :, :, :] .=
                    getfield(rays, field)[1:local_count, i0:i1, j0:j1, k1]
            end
            MPI.Send(rays_up, comm; dest = up)
        end
    end

    return
end

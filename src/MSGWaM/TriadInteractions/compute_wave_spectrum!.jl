"""
```julia
compute_wave_spectrum!(state::State)
```

Compute the gravity-wave spectrum needed for the computation of the scattering integral by dispatching to a WKB-mode-specific method.

```julia
compute_wave_spectrum!(state::State, wkb_mode::Union{MultiColumn, SingleColumn})
```

Compute the gravity-wave spectrum needed for the computation of the scattering integral in multi-column mode.

This method computes the sums


Therein, ``\\left(\\lambda, \\mu, \\nu\\right)`` are index shifts to ray volumes that are at least partially within the grid cell at ``\\left(i, j, k\\right)``, ``F_{r, i + \\lambda, j + \\mu, k + \\nu}`` are the corresponding ray volume fractions and ``\\left(u_\\mathrm{w}, v_\\mathrm{w}, w_\\mathrm{w}, \\theta_\\mathrm{w}\\right)_{r, i + \\lambda, j + \\mu, k + \\nu}`` are the wave amplitudes of the wind and the potential temperature. The computation is based on the relations

```math
\\begin{align*}
    \\overline{\\rho} u_{\\mathrm{w}, r} u_{\\mathrm{w}, r}^* & = \\left(k_r \\widehat{c}_{\\mathrm{g} x, r} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_r \\widehat{c}_{\\mathrm{g} x, r} + l_r \\widehat{c}_{\\mathrm{g} y, r}}{1 - \\left(\\widehat{\\omega}_r / f\\right)^2}\\right) \\mathcal{A}_r,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = l_r \\widehat{c}_{\\mathrm{g} x, r} \\mathcal{A}_r,\\\\
    \\overline{\\rho} u_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = \\frac{k_r \\widehat{c}_{\\mathrm{g} z, r}}{1 - \\left(f / \\widehat{\\omega}_r\\right)^2} \\mathcal{A}_r,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = \\left(l_r \\widehat{c}_{\\mathrm{g} y, r} - \\mathrm{sgn} \\left(\\left|f\\right|\\right) \\frac{k_r \\widehat{c}_{\\mathrm{g} x, r} + l_r \\widehat{c}_{\\mathrm{g} y, r}}{1 - \\left(\\widehat{\\omega}_r / f\\right)^2}\\right) \\mathcal{A}_r,\\\\
    \\overline{\\rho} v_{\\mathrm{w}, r} w_{\\mathrm{w}, r}^* & = \\frac{l_r \\widehat{c}_{\\mathrm{g} z, r}}{1 - \\left(f / \\widehat{\\omega}_r\\right)^2} \\mathcal{A}_r,\\\\
    \\theta_{\\mathrm{w}, r} u_{\\mathrm{w}, r}^* & = \\frac{f \\overline{\\theta}}{g \\overline{\\rho}} \\frac{l_r m_r N_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r,\\\\
    \\theta_{\\mathrm{w}, r} v_{\\mathrm{w}, r}^* & = - \\frac{f \\overline{\\theta}}{g \\overline{\\rho}} \\frac{k_r m_r N_r^2}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{A}_r,
\\end{align*}
```

where ``N_r^2`` is the squared buoyancy frequency interpolated to the ray-volume position. The components of the intrinsic group velocity are given by

```math
\\begin{align*}
    \\widehat{c}_{\\mathrm{g} x, r} & = \\frac{k_r \\left(N_r^2 - \\widehat{\\omega}_r^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} y, r} & = \\frac{l_r \\left(N_r^2 - \\widehat{\\omega}_r^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2},\\\\
    \\widehat{c}_{\\mathrm{g} z, r} & = - \\frac{m_r \\left(\\widehat{\\omega}_r^2 - f^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2}.
\\end{align*}
```

# Arguments

  - `state::State`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.MSGWaM.TriadInteractions.compute_spectral_cell_indices`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
"""
function compute_wave_spectrum! end

function compute_wave_spectrum!(state::State)
    (; wkb_mode) = state.namelists.wkb
    (; triad_int) = state.namelists.triad
    if triad_int
        compute_wave_spectrum!(state, wkb_mode, triad_int)
    end
    return
end

function compute_wave_spectrum!(state::State, wkb_mode::Union{MultiColumn, SingleColumn}, triad_int::Bool)
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branch) = state.namelists.wkb
    (; tref, g_ndim) = state.constants
    (; i0, i1, j0, j1, k0, k1) = domain
    (; kp, m, dx, dy, dz, x, y, zctilde, jac) = grid
    (; nray, rays, spec_tend) = state.wkb

    (ukp, lkp) = half_logwidth(kp)
    (um, lm) = half_logwidth(m)
    

    # Set Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(TriadTendencies)
        getfield(spec_tend, field) .= 0.0
    end



    @ivy for k in (k0 - 1):(k1 + 1),
        j in (j0 - 1):(j1 + 1),
        i in (i0 - 1):(i1 + 1)

        for r in 1:nray[i, j, k]
            if rays.dens[r, i, j, k] == 0
                continue
            end

            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            dxr = rays.dxray[r, i, j, k]
            dyr = rays.dyray[r, i, j, k]
            dzr = rays.dzray[r, i, j, k]

            kr = rays.k[r, i, j, k]
            lr = rays.l[r, i, j, k]
            mr = abs(rays.m[r, i, j, k])

            dkr = rays.dkray[r, i, j, k]
            dlr = rays.dlray[r, i, j, k]
            dmr = abs(rays.dmray[r, i, j, k])

            kpr = sqrt(kr^2 + lr^2)
            dkpr = sqrt(dkr^2 + dlr^2)

            (imin, imax, jmin, jmax) =
                compute_horizontal_cell_indices(state, xr, yr, dxr, dyr)

            (kpmin, kpmax, mmin, mmax) = compute_spectral_cell_indices(state, kpr, mr, dkpr, dmr)
            
            

            for iray in imin:imax
                if x_size > 1
                    dxi = (
                        min(xr + dxr / 2, x[iray] + dx / 2) -
                        max(xr - dxr / 2, x[iray] - dx / 2)
                    )

                    fcpspx =  dxi / dx
                else
                    fcpspx = 1.0
                end

                for jray in jmin:jmax
                    if y_size > 1
                        dyi = (
                            min(yr + dyr / 2, y[jray] + dy / 2) -
                            max(yr - dyr / 2, y[jray] - dy / 2)
                        )

                        fcpspy =  dyi / dy
                    else
                        fcpspy = 1.0
                    end

                    kmin = get_next_half_level(
                        iray,
                        jray,
                        zr - dzr / 2,
                        state;
                        dkd = 1,
                    )
                    kmax = get_next_half_level(
                        iray,
                        jray,
                        zr + dzr / 2,
                        state;
                        dkd = 1,
                    )

                    for kray in kmin:kmax
                        dzi =
                            min((zr + dzr / 2), zctilde[iray, jray, kray]) -
                            max((zr - dzr / 2), zctilde[iray, jray, kray - 1])

                        fcpspz =  dzi / jac[iray, jray, kray] / dz
                        
                        for kpray in kpmin:kpmax
                             dkpi = 
                                min(kr + dkr / 2, kp[kpray] + ukp[kpray]) -
                                max(kr - dkr / 2, kp[kpray] - lkp[kpray])
                                      
                             fcpspkp =  dkpi / dkpr

                             for mray in mmin:mmax
                                dmi = 
                                    min(mr + dmr / 2, m[mray] + um[mray]) -
                                    max(mr - dmr / 2, m[mray] - lm[mray])
                                fcpspm = dmi / dmr
                                
                                wadr = fcpspx * fcpspy * fcpspz * fcpspkp * fcpspm * rays.dens[r, i, j, k]
                                spec_tend.wavespectrum[iray, jray, kray, kpray, mray] += wadr

                             end

                        end

                    end

                end

            end

        end
        

    end
    

end


                        
            
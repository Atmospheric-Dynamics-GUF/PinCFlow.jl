function launch_new_ray_vol! end
#=
function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad2D)

     (; nray, rays, nray_wrk) = state.wkb
     (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid
     (;lref ) = state.constants

     
     ray_index = nray[i, j, k] + 1
    
     @ivy if ray_index > nray_wrk
            error(
                "Error in lanching new ray volumes!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
    end


     rays.x[ray_index, i, j, k] = x[i]
     rays.y[ray_index, i, j, k] = y[j]
     rays.z[ray_index, i, j, k] = zc[i, j, k]
     rays.k[ray_index, i, j, k] =  kpr ##right now it is launching the rays with positive wave numbers only.
     rays.l[ray_index, i, j, k] = 0.0 ##need to include the negative wave numbers also
     rays.m[ray_index, i, j, k] = mr
     rays.dxray[ray_index, i, j, k] = dx
     rays.dyray[ray_index, i, j, k] = dy
     rays.dzray[ray_index, i, j, k] = dz
     rays.dkray[ray_index, i, j, k] = dkpr
     rays.dlray[ray_index, i, j, k] = 0.0
     rays.dmray[ray_index, i, j, k] = dmr
     rays.dens[ray_index, i, j, k] = was
     nray[i, j, k] += 1

    #println("\n with the non dimensional ray volume properties ", 
    #            (kpr, mr, dkpr, dmr, was))
end
=#
function launch_new_ray_vol!(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    kpr::AbstractFloat,
    mr::AbstractFloat,
    dkpr::AbstractFloat,
    dmr::AbstractFloat,
    was::AbstractFloat,
    triad_mode::Triad2D,
)

    (; nray, rays, nray_wrk, cgx_max, cgz_max) = state.wkb
    (; dxray, dyray, dzray, dkray, dlray, dmray, ddxray, ddyray, ddzray) = state.wkb.increments
    (; dx, dy, dz, x, y, zc, jac) = state.grid
    (; tref) = state.constants
    (; x_size) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; nrx, nrz, branch, wkb_mode) = state.namelists.wkb

    fc = coriolis_frequency * tref

    # For x_size == 1 there is no physical subdivision in x.
    nrx_launch = x_size == 1 ? 1 : nrx
    nray_new = nrx_launch * nrz

    # Make sure sufficient ray-volume storage is available before
    # launching any new ray volumes.
    if nray[i, j, k] + nray_new > nray_wrk
        error(
            "Error in launching new ray volumes!: required nray",
            [i, j, k],
            " = ",
            nray[i, j, k] + nray_new,
            " > nray_wrk = ",
            nray_wrk,
        )
    end

    ray_index = nray[i, j, k] + 1

    for ix in 1:nrx_launch, kz in 1:nrz

        #----------------------------------------------------------
        # Physical position and extent
        #----------------------------------------------------------

        if x_size == 1
            rays.x[ray_index, i, j, k] = x[i]
            rays.dxray[ray_index, i, j, k] = dx
        else
            rays.x[ray_index, i, j, k] =
                x[i] - 0.5 * dx + (ix - 0.5) * dx / nrx_launch

            rays.dxray[ray_index, i, j, k] = dx / nrx_launch
        end

        rays.y[ray_index, i, j, k] = y[j]

        rays.z[ray_index, i, j, k] =
            zc[i, j, k] - 0.5 * jac[i, j, k] * dz +
            (kz - 0.5) * jac[i, j, k] * dz / nrz

        rays.dyray[ray_index, i, j, k] = dy
        rays.dzray[ray_index, i, j, k] = jac[i, j, k] * dz / nrz

        #----------------------------------------------------------
        # Spectral position and extent
        #----------------------------------------------------------

        rays.k[ray_index, i, j, k] = kpr
        rays.l[ray_index, i, j, k] = 0.0
        rays.m[ray_index, i, j, k] = mr

        # For x_size == 1, kp is a discrete mode label and there
        # is no physical spectral extent in the k-direction.
        rays.dkray[ray_index, i, j, k] = x_size == 1 ? 0.0 : dkpr
        rays.dlray[ray_index, i, j, k] = 0.0
        rays.dmray[ray_index, i, j, k] = dmr

        # Wave-action density.
        rays.dens[ray_index, i, j, k] = was

        #----------------------------------------------------------
        # Initialize Runge-Kutta increments of the new ray volume
        #----------------------------------------------------------

        dxray[ray_index, i, j, k] = 0.0
        dyray[ray_index, i, j, k] = 0.0
        dzray[ray_index, i, j, k] = 0.0

        dkray[ray_index, i, j, k] = 0.0
        dlray[ray_index, i, j, k] = 0.0
        dmray[ray_index, i, j, k] = 0.0

        ddxray[ray_index, i, j, k] = 0.0
        ddyray[ray_index, i, j, k] = 0.0
        ddzray[ray_index, i, j, k] = 0.0

        #----------------------------------------------------------
        # Update maximum group velocities used for the WKB CFL
        #----------------------------------------------------------

        xr = rays.x[ray_index, i, j, k]
        yr = rays.y[ray_index, i, j, k]
        zr = rays.z[ray_index, i, j, k]

        n2r = interpolate_stratification(zr, state, N2())

        if n2r < 0.0
            error("Error in launch_new_ray_vol!: interpolated stratification is negative.")
        end

        wnrk = rays.k[ray_index, i, j, k]
        wnrl = rays.l[ray_index, i, j, k]
        wnrm = rays.m[ray_index, i, j, k]

        wnrh = sqrt(wnrk^2 + wnrl^2)

        if wnrh <= 0.0
            error("Error in launch_new_ray_vol!: horizontal wavenumber is non-positive.")
        end

        # Use the same non-hydrostatic dispersion relation as the
        # MS-GWaM ray-propagation routine.
        omir =
            branch *
            sqrt(n2r * wnrh^2 + fc^2 * wnrm^2) /
            sqrt(wnrh^2 + wnrm^2)

        # Horizontal group velocity contributes to the WKB CFL
        # only when horizontal ray propagation is active.
        if x_size > 1 && wkb_mode != SingleColumn()
            uxr = interpolate_mean_flow(xr, yr, zr, state, U())

            cgirx =
                wnrk * (n2r - omir^2) /
                (omir * (wnrh^2 + wnrm^2))

            cgx_max[] = max(cgx_max[], abs(uxr + cgirx))
        end

        # Vertical group velocity. This follows the convention used
        # in propagate_rays!, where no resolved vertical mean-flow
        # velocity is added to c_gz.
        cgirz =
            -wnrm * (omir^2 - fc^2) /
            (omir * (wnrh^2 + wnrm^2))

        cgz_max[i, j, k] = max(cgz_max[i, j, k], abs(cgirz))

        # Ray volume has now been completely initialized.
        ray_index += 1
        nray[i, j, k] += 1
    end

    return nothing
end


function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad3DIso)

     (; nray, rays, nray_wrk) = state.wkb
     (; dx, dy, dz, x, y, zc ,zctilde, jac) = state.grid

     
     ray_index = nray[i, j, k] + 1
    
     @ivy if ray_index > nray_wrk
            error(
                "Error in lanching new ray volumes!: nray",
                [i, j, k],
                " > nray_wrk = ",
                nray_wrk,
            )
    end


     rays.x[ray_index, i, j, k] = x[i]
     rays.y[ray_index, i, j, k] = y[j]
     rays.z[ray_index, i, j, k] = zc[i, j, k]
     rays.k[ray_index, i, j, k] = sqrt(2) * kpr ##right now it is launching the rays with positive wave numbers only.
     rays.l[ray_index, i, j, k] = sqrt(2) * kpr ##need to include the negative wave numbers also
     rays.m[ray_index, i, j, k] = mr
     rays.dxray[ray_index, i, j, k] = dx
     rays.dyray[ray_index, i, j, k] = dy
     rays.dzray[ray_index, i, j, k] = dz
     rays.dkray[ray_index, i, j, k] = sqrt(2) * dkpr
     rays.dlray[ray_index, i, j, k] = sqrt(2) * dkpr
     rays.dmray[ray_index, i, j, k] = dmr
     rays.dens[ray_index, i, j, k] = was
     nray[i, j, k] += 1


end


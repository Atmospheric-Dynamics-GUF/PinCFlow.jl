function launch_new_ray_vol! end

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

function launch_new_ray_vol!(state::State, i::Integer, j::Integer, k::Integer, kpr::AbstractFloat, mr::AbstractFloat,
     dkpr::AbstractFloat, dmr::AbstractFloat, was::AbstractFloat, triad_mode::Triad2D)

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


end
function check_ray_volumes(state::State)
  (; sizex, sizey) = state.namelists.domain
  (; io, jo, i0, i1, j0, j1, k0, k1) = state.state.domain
  (; x, y, ztildetfc, dx, dy)
  (; nray, rays) = state.wkb

  # Loop over ray volumes.
  for kz in k0:k1, jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
    if nray[ix, jy, kz] > 0
      for iray in 1:nray[ix, jy, kz]

        # Check zonal position.
        if sizex > 1
          xr = rays.x[iray, ix, jy, kz]

          if xr < x[io + ix] - dx / 2
            println("Error in shift_ray_volumes!:")
            println(
              "xr = ",
              xr,
              "< x(io + ix) - dx / 2 = ",
              x[io + ix] - dx / 2,
            )
            println("io = ", io)
            println(
              "(iRay, ix, jy, kz) = (",
              iray,
              ", ",
              ix,
              ", ",
              jy,
              ", ",
              kz,
              ")",
            )
            exit()
          end

          if xr > x[io + ix] + dx / 2
            println("Error in shift_ray_volumes!:")
            println(
              "xr = ",
              xr,
              "> x(io + ix) + dx / 2 = ",
              x[io + ix] + dx / 2,
            )
            println("io = ", io)
            println(
              "(iRay, ix, jy, kz) = (",
              iray,
              ", ",
              ix,
              ", ",
              jy,
              ", ",
              kz,
              ")",
            )
            exit()
          end
        end

        # Check meridional position.
        if sizey > 1
          yr = rays.y[iray, ix, jy, kz]

          if yr < y[jo + jy] - dy / 2
            println("Error in shift_ray_volumes!:")
            println(
              "yr = ",
              yr,
              "< y(jo + jy) - dy / 2 = ",
              y[jo + jy] - dy / 2,
            )
            println("jo = ", jo)
            println(
              "(iRay, ix, jy, kz) = (",
              iray,
              ", ",
              ix,
              ", ",
              jy,
              ", ",
              kz,
              ")",
            )
            exit()
          end

          if yr > y[jo + jy] + dy / 2
            println("Error in shift_ray_volumes!:")
            println(
              "yr = ",
              yr,
              "> y(jo + jy) + dy / 2 = ",
              y[jo + jy] + dy / 2,
            )
            println("jo = ", jo)
            println(
              "(iRay, ix, jy, kz) = (",
              iray,
              ", ",
              ix,
              ", ",
              jy,
              ", ",
              kz,
              ")",
            )
            exit()
          end
        end

        # Check vertical position.
        zr = rays.z[iray, ix, jy, kz]

        if zr < ztildetfc[ix, jy, kz - 1]
          println("Error in shift_ray_volumes!:")
          println(
            "zr =",
            zr,
            "< ztildetfc(ix, jy, kz - 1) = ",
            ztildetfc[ix, jy, kz - 1],
          )
          println(
            "(iray, ix, jy, kz) = (",
            iray,
            ", ",
            ix,
            ", ",
            jy,
            ", ",
            kz,
            ")",
          )
          exit()
        end

        if zr > ztildetfc[ix, jy, kz]
          println("Error in shift_ray_volumes!:")
          println(
            "zr =",
            zr,
            "> ztildetfc(ix, jy, kz) = ",
            ztildetfc[ix, jy, kz],
          )
          println(
            "(iray, ix, jy, kz) = (",
            iray,
            ", ",
            ix,
            ", ",
            jy,
            ", ",
            kz,
            ")",
          )
          exit()
        end
      end
    end
  end
  return
end
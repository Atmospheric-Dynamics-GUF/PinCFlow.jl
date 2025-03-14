struct Grid{
  A <: OffsetVector{<:AbstractFloat},
  B <: AbstractFloat,
  C <: OffsetMatrix{<:AbstractFloat},
  D <: OffsetArray{<:AbstractFloat, 3},
  E <: OffsetArray{<:AbstractFloat, 5},
}

  # Scaled domain.
  lx::A
  ly::A
  lz::A

  # Grid spacings.
  dx::B
  dy::B
  dz::B

  # Coordinates.
  x::A
  y::A
  z::A

  # Stretched vertical grid.
  zs::A
  ztildes::A

  # Topography.
  topography_surface::C

  # Jacobian and metric tensor.
  jac::D
  met::E

  # Vertical layers.
  ztfc::D
  ztildetfc::D
end

function Grid(namelists::Namelists, constants::Constants, domain::Domain)

  # Get parameters.
  (; sizex, sizey, sizez, lx_dim, ly_dim, lz_dim, nbx, nby, nbz) =
    namelists.domain
  (;
    mountainheight_dim,
    mountainwidth_dim,
    mountain_case,
    range_factor,
    spectral_modes,
    envelope_reduction,
    stretch_exponent,
  ) = namelists.grid
  (; is, js, nx, ny, nz) = domain
  (; lref) = constants

  # Non-dimensionalize domain boundaries.
  lx = OffsetVector(lx_dim, 0:1) ./ lref
  ly = OffsetVector(ly_dim, 0:1) ./ lref
  lz = OffsetVector(lz_dim, 0:1) ./ lref

  # Compute grid spacings.
  dx = (lx[1] - lx[0]) / sizex
  dy = (ly[1] - ly[0]) / sizey
  dz = (lz[1] - lz[0]) / sizez

  # Compute x-coordinate.
  x = OffsetVector(zeros(sizex + 1 + 2 * nbx), (-nbx):(sizex + nbx))
  for ix in (-nbx):(sizex + nbx)
    x[ix] = lx[0] + (ix - 1) * dx + dx / 2
  end

  # Compute y-coordinate.
  y = OffsetVector(zeros(sizey + 1 + 2 * nby), (-nby):(sizey + nby))
  for jy in (-nby):(sizey + nby)
    y[jy] = ly[0] + (jy - 1) * dy + dy / 2
  end

  # Compute z-coordinate.
  z = OffsetVector(zeros(sizez + 1 + 2 * nbz), (-nbz):(sizez + nbz))
  for kz in (-nbz):(sizez + nbz)
    z[kz] = lz[0] + (kz - 1) * dz + dz / 2
  end

  # Initialize the stretched vertical grid.
  (ztildes, zs) =
    (OffsetVector(zeros(nz + 2 * nbz + 1), (-nbz):(nz + nbz)) for i in 1:2)

  # Compute the stretched vertical grid.
  for kz in (-nbz):(nz + nbz)
    level = z[kz] + 0.5 * dz
    if level < 0
      ztildes[kz] = -lz[1] * (-level / lz[1])^stretch_exponent
    elseif level > lz[1]
      ztildes[kz] =
        2 * lz[1] - lz[1] * ((2 * lz[1] - level) / lz[1])^stretch_exponent
    else
      ztildes[kz] = lz[1] * (level / lz[1])^stretch_exponent
    end
  end
  for kz in (-nbz + 1):(nz + nbz)
    zs[kz] = 0.5 * (ztildes[kz] + ztildes[kz - 1])
  end
  zs[-nbz] = ztildes[-nbz] - 0.5 * (ztildes[nbz + 1] - ztildes[nbz])

  if lz[0] != 0.0
    println("Error in setup_topography: lz[0] must be zero for TFC!")
    exit()
  end

  mountainheight = mountainheight_dim / lref
  mountainwidth = mountainwidth_dim / lref
  mountainwavenumber = pi / mountainwidth

  x_center = 0.5 * (lx[1] + lx[0])
  y_center = 0.5 * (ly[1] + ly[0])

  ix0 = is + nbx - 1
  jy0 = js + nby - 1

  topography_surface = OffsetMatrix(
    zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1)),
    (-nbx):(nx + nbx),
    (-nby):(ny + nby),
  )

  if mountain_case != 0
    for jy in 1:(ny)
      for ix in 1:(nx)
        if mountain_case == 1
          # 2D cosine mountains
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            (1.0 + cos(mountainwavenumber * (x[ix + ix0] - x_center)))

        elseif mountain_case == 2
          # 3D cosine mountains
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            (
              1.0 + cos(
                mountainwavenumber * sqrt(
                  (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                ),
              )
            )

        elseif mountain_case == 3
          # 2D isolated mountain
          topography_surface[ix, jy] =
            mountainheight /
            (1.0 + (x[ix + ix0] - x_center)^2.0 / mountainwidth^2.0)

        elseif mountain_case == 4
          # 3D isolated mountain
          topography_surface[ix, jy] =
            mountainheight / (
              1.0 +
              ((x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0) /
              mountainwidth^2.0
            )

        elseif mountain_case == 5
          # 2D cosine envelope and even background
          if abs(x[ix + ix0] - x_center) <= mountainwidth * range_factor
            topography_surface[ix, jy] =
              0.5 *
              mountainheight *
              (
                1.0 +
                0.5 *
                (
                  1.0 + cos(
                    mountainwavenumber / range_factor *
                    (x[ix + ix0] - x_center),
                  )
                ) *
                cos(mountainwavenumber * (x[ix + ix0] - x_center))
              )
          else
            topography_surface[ix, jy] = 0.5 * mountainheight
          end

        elseif mountain_case == 6
          # 3D cosine envelope and even background
          if (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0 <=
             (mountainwidth * range_factor)^2.0
            topography_surface[ix, jy] =
              0.5 *
              mountainheight *
              (
                1.0 +
                0.5 *
                (
                  1.0 + cos(
                    mountainwavenumber / range_factor * sqrt(
                      (x[ix + ix0] - x_center)^2.0 +
                      (y[jy + jy0] - y_center)^2.0,
                    ),
                  )
                ) *
                cos(
                  mountainwavenumber * sqrt(
                    (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                  ),
                )
              )
          else
            topography_surface[ix, jy] = 0.5 * mountainheight
          end

        elseif mountain_case == 7
          # 2D Gaussian envelope and even background
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            (
              1.0 +
              exp(
                -(x[ix + ix0] - x_center)^2.0 /
                (mountainwidth * range_factor)^2.0,
              ) * cos(mountainwavenumber * (x[ix + ix0] - x_center))
            )

        elseif mountain_case == 8
          # 3D Gaussian envelope and even background
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            (
              1.0 +
              exp(
                -((x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0) /
                (mountainwidth * range_factor)^2.0,
              ) * cos(
                mountainwavenumber * sqrt(
                  (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                ),
              )
            )

        elseif mountain_case == 9
          # 2D cosine envelope and cosine background
          if abs(x[ix + ix0] - x_center) <= mountainwidth * range_factor
            topography_surface[ix, jy] =
              0.25 *
              mountainheight *
              (
                1.0 + cos(
                  mountainwavenumber / range_factor * (x[ix + ix0] - x_center),
                )
              ) *
              (1.0 + cos(mountainwavenumber * (x[ix + ix0] - x_center)))
          end

        elseif mountain_case == 10
          # 3D cosine envelope and cosine background
          if (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0 <=
             (mountainwidth * range_factor)^2.0
            topography_surface[ix, jy] =
              0.25 *
              mountainheight *
              (
                1.0 + cos(
                  mountainwavenumber / range_factor * sqrt(
                    (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                  ),
                )
              ) *
              (
                1.0 + cos(
                  mountainwavenumber * sqrt(
                    (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                  ),
                )
              )
          end

        elseif mountain_case == 11
          # 2D Gaussian envelope and Gaussian background
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            exp(
              -(x[ix + ix0] - x_center)^2.0 /
              (mountainwidth * range_factor)^2.0,
            ) *
            (1.0 + cos(mountainwavenumber * (x[ix + ix0] - x_center)))

        elseif mountain_case == 12
          # 3D Gaussian envelope and Gaussian background
          topography_surface[ix, jy] =
            0.5 *
            mountainheight *
            exp(
              -((x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0) /
              (mountainwidth * range_factor)^2.0,
            ) *
            (
              1.0 + cos(
                mountainwavenumber * sqrt(
                  (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                ),
              )
            )

        elseif mountain_case == 13
          # 3D WKB topography
          if (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0 <=
             (mountainwidth * range_factor)^2.0
            topography_surface[ix, jy] =
              0.25 *
              mountainheight *
              (
                1.0 + cos(
                  mountainwavenumber / range_factor * sqrt(
                    (x[ix + ix0] - x_center)^2.0 + (y[jy + jy0] - y_center)^2.0,
                  ),
                )
              ) *
              (1.0 + envelope_reduction)
            for iwm in 0:(spectral_modes - 1)
              kk = mountainwavenumber * cos(pi / spectral_modes * iwm)
              ll = mountainwavenumber * sin(pi / spectral_modes * iwm)
              topography_surface[ix, jy] =
                topography_surface[ix, jy] +
                0.25 *
                mountainheight *
                (
                  1.0 + cos(
                    mountainwavenumber / range_factor * sqrt(
                      (x[ix + ix0] - x_center)^2.0 +
                      (y[jy + jy0] - y_center)^2.0,
                    ),
                  )
                ) *
                cos(
                  kk * (x[ix + ix0] - x_center) + ll * (y[jy + jy0] - y_center),
                ) / spectral_modes * (1.0 - envelope_reduction)
            end
          end
        else
          println("Error in setup_topography: Unknown mountain case!")
          exit()
        end
      end
    end
  else
    println("Error in setup_topography: Mountain case 0 is not ready yet!")
  end

  # Set halos of topography surface.
  set_zonal_boundaries_of_field!(topography_surface, namelists, domain)
  set_meridional_boundaries_of_field!(topography_surface, namelists, domain)

  # Initialize Jacobian and metric tensor.
  jac = OffsetArray(
    zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
    (-nbx):(nx + nbx),
    (-nby):(ny + nby),
    (-nbz):(nz + nbz),
  )
  met = OffsetArray(
    zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1, 3, 3)),
    (-nbx):(nx + nbx),
    (-nby):(ny + nby),
    (-nbz):(nz + nbz),
    1:3,
    1:3,
  )

  # Compute the Jacobian.
  for kz in (-nbz + 1):(nz + nbz)
    jac[:, :, kz] =
      (lz[1] .- topography_surface) ./ lz[1] .*
      (ztildes[kz] .- ztildes[kz - 1]) ./ dz
  end
  jac[:, :, -nbz] = jac[:, :, nbz + 1]

  # Compute the metric tensor.
  met[:, :, :, 1, 1] .= 1.0
  met[:, :, :, 1, 2] .= 0.0
  for kz in (-nbz + 1):(nz + nbz)
    for jy in 1:(ny)
      for ix in 1:(nx)
        met[ix, jy, kz, 1, 3] =
          (topography_surface[ix + 1, jy] - topography_surface[ix - 1, jy]) /
          (2.0 * dx) * (zs[kz] - lz[1]) / (lz[1] - topography_surface[ix, jy]) *
          dz / (ztildes[kz] - ztildes[kz - 1])
      end
    end
  end
  set_zonal_boundaries_of_field!(met[:, :, :, 1, 3], namelists, domain)
  set_meridional_boundaries_of_field!(met[:, :, :, 1, 3], namelists, domain)
  met[:, :, -nbz, 1, 3] =
    met[:, :, nbz + 1, 1, 3] .* (zs[-nbz] .- lz[1]) ./ (zs[nbz + 1] .- lz[1])
  met[:, :, :, 2, 1] .= 0.0
  met[:, :, :, 2, 2] .= 1.0
  for kz in (-nbz + 1):(nz + nbz)
    for jy in 1:(ny)
      for ix in 1:(nx)
        met[ix, jy, kz, 2, 3] =
          (topography_surface[ix, jy + 1] - topography_surface[ix, jy - 1]) /
          (2.0 * dy) * (zs[kz] - lz[1]) / (lz[1] - topography_surface[ix, jy]) *
          dz / (ztildes[kz] - ztildes[kz - 1])
      end
    end
  end
  set_zonal_boundaries_of_field!(met[:, :, :, 2, 3], namelists, domain)
  set_meridional_boundaries_of_field!(met[:, :, :, 2, 3], namelists, domain)
  met[:, :, -nbz, 2, 3] =
    met[:, :, nbz + 1, 2, 3] .* (zs[-nbz] .- lz[1]) ./ (zs[nbz + 1] .- lz[1])
  met[:, :, :, 3, 1] = met[:, :, :, 1, 3]
  met[:, :, :, 3, 2] = met[:, :, :, 2, 3]
  for kz in (-nbz + 1):(nz + nbz)
    for jy in 1:(ny)
      for ix in 1:(nx)
        met[ix, jy, kz, 3, 3] =
          (
            (lz[1] / (lz[1] - topography_surface[ix, jy]))^2.0 +
            ((zs[kz] - lz[1]) / (lz[1] - topography_surface[ix, jy]))^2.0 * (
              (
                (
                  topography_surface[ix + 1, jy] -
                  topography_surface[ix - 1, jy]
                ) / (2.0 * dx)
              )^2.0 +
              (
                (
                  topography_surface[ix, jy + 1] -
                  topography_surface[ix, jy - 1]
                ) / (2.0 * dy)
              )^2.0
            )
          ) * (dz / (ztildes[kz] - ztildes[kz - 1]))^2.0
      end
    end
  end
  set_zonal_boundaries_of_field!(met[:, :, :, 3, 3], namelists, domain)
  set_meridional_boundaries_of_field!(met[:, :, :, 3, 3], namelists, domain)
  for jy in 1:(ny)
    for ix in 1:(nx)
      met[ix, jy, -nbz, 3, 3] =
        (
          (lz[1] / (lz[1] - topography_surface[ix, jy]))^2.0 +
          ((zs[-nbz] - lz[1]) / (lz[1] - topography_surface[ix, jy]))^2.0 * (
            (
              (
                topography_surface[ix + 1, jy] - topography_surface[ix - 1, jy]
              ) / (2.0 * dx)
            )^2.0 +
            (
              (
                topography_surface[ix, jy + 1] - topography_surface[ix, jy - 1]
              ) / (2.0 * dy)
            )^2.0
          )
        ) * (dz / (ztildes[nbz + 1] - ztildes[nbz]))^2.0
    end
  end
  set_zonal_boundaries_of_field!(met[:, :, -nbz, 3, 3], namelists, domain)
  set_meridional_boundaries_of_field!(met[:, :, -nbz, 3, 3], namelists, domain)

  # Initialize the physical layers.
  (ztildetfc, ztfc) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:2
  )

  # Compute the physical layers.
  for kz in (-nbz):(nz + nbz)
    ztildetfc[:, :, kz] =
      (lz[1] .- topography_surface) ./ lz[1] .* ztildes[kz] .+
      topography_surface
    ztfc[:, :, kz] =
      (lz[1] .- topography_surface) ./ lz[1] .* zs[kz] .+ topography_surface
  end

  # Return a Grid instance.
  return Grid(
    lx,
    ly,
    lz,
    dx,
    dy,
    dz,
    x,
    y,
    z,
    zs,
    ztildes,
    topography_surface,
    jac,
    met,
    ztfc,
    ztildetfc,
  )
end

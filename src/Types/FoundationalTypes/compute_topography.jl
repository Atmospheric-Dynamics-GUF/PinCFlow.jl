"""
```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::WKBMountainWave,
)
```

Compute topography for WKB mountain wave test cases with spectral decomposition.

# Arguments

  - `namelists::Namelists`: Configuration parameters including domain, grid, and WKB settings
  - `constants::Constants`: Physical constants including reference length scale
  - `domain::Domain`: Computational domain specification
  - `x::AbstractVector{<:AbstractFloat}`: X-coordinate grid points
  - `y::AbstractVector{<:AbstractFloat}`: Y-coordinate grid points
  - `testcase::WKBMountainWave`: WKB mountain wave test case specification

# Returns

  - `Tuple`: (topography_surface, topography_spectrum, k_spectrum, l_spectrum)

      + `topography_surface`: 2D array of surface height values
      + `topography_spectrum`: 3D array of spectral amplitude components
      + `k_spectrum`: 3D array of zonal wavenumber components
      + `l_spectrum`: 3D array of meridional wavenumber components

# Mountain Cases (WKB-specific)

  - **Case 1**: 2D cosine mountains - Simple periodic structure in x-direction
  - **Case 5**: 2D cosine envelope with even background - Modulated cosine wave
  - **Case 7**: 2D Gaussian envelope with even background - Gaussian-modulated structure
  - **Case 9**: 2D cosine envelope with cosine background - Double cosine modulation
  - **Case 11**: 2D Gaussian envelope with Gaussian background - Full Gaussian structure
  - **Case 13**: 3D WKB topography - Complex 3D structure with multiple spectral modes

# Notes

  - For Case 13, spectral modes are distributed uniformly in wavenumber space
  - The function validates that `nwm` (number of wave modes) is sufficient for the chosen case
  - All coordinates are normalized by the reference length scale
"""
function compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::WKBMountainWave,
)
    (; lx_dim, ly_dim) = namelists.domain
    (; testcase) = namelists.setting
    (; nwm) = namelists.wkb
    (;
        mountainheight_dim,
        mountainwidth_dim,
        mountain_case,
        height_factor,
        width_factor,
        spectral_modes,
    ) = namelists.grid
    (; nxx, nyy, io, jo, i0, i1, j0, j1) = domain
    (; lref) = constants

    if nwm < 1 || (mountain_case == 13 && nwm < spectral_modes)
        error("Error in compute_topography: nwm is too small!")
    end

    mountainheight = mountainheight_dim / lref
    mountainwidth = mountainwidth_dim / lref
    mountainwavenumber = pi / mountainwidth

    x_center = 0.5 * (lx_dim[2] + lx_dim[1]) / lref
    y_center = 0.5 * (ly_dim[2] + ly_dim[1]) / lref

    topography_surface = zeros(nxx, nyy)
    topography_spectrum = zeros(nwm, nxx, nyy)
    k_spectrum = zeros(nwm, nxx, nyy)
    l_spectrum = zeros(nwm, nxx, nyy)

    for jy in j0:j1, ix in i0:i1
        if mountain_case == 1
            # 2D cosine mountains
            topography_surface[ix, jy] = 0.5 * mountainheight
            topography_spectrum[1, ix, jy] = 0.5 * mountainheight
            k_spectrum[1, ix, jy] = mountainwavenumber

        elseif mountain_case == 5
            # 2D cosine envelope and even background
            if abs(x[io + ix] - x_center) <= mountainwidth * width_factor
                k_spectrum[1, ix, jy] = mountainwavenumber
                topography_spectrum[1, ix, jy] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            (x[io + ix] - x_center),
                        )
                    )
            end
            topography_surface[ix, jy] = 0.5 * mountainheight

        elseif mountain_case == 7
            # 2D Gaussian envelope and even background
            k_spectrum[1, ix, jy] = mountainwavenumber
            topography_spectrum[1, ix, jy] =
                0.5 *
                mountainheight *
                exp(
                    -(x[io + ix] - x_center)^2.0 /
                    (mountainwidth * width_factor)^2.0,
                )
            topography_surface[ix, jy] = 0.5 * mountainheight

        elseif mountain_case == 9
            # 2D cosine envelope and cosine background
            if abs(x[io + ix] - x_center) <= mountainwidth * width_factor
                k_spectrum[1, ix, jy] = mountainwavenumber
                topography_spectrum[1, ix, jy] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            (x[io + ix] - x_center),
                        )
                    )
                topography_surface[ix, jy] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            (x[io + ix] - x_center),
                        )
                    )
            end

        elseif mountain_case == 11
            # 2D Gaussian envelope and Gaussian background
            k_spectrum[1, ix, jy] = mountainwavenumber
            topography_spectrum[1, ix, jy] =
                0.5 *
                mountainheight *
                exp(
                    -(x[io + ix] - x_center)^2.0 /
                    (mountainwidth * width_factor)^2.0,
                )
            topography_surface[ix, jy] =
                0.5 *
                mountainheight *
                exp(
                    -(x[io + ix] - x_center)^2.0 /
                    (mountainwidth * width_factor)^2.0,
                )

        elseif mountain_case == 13
            # 3D WKB topography
            if (x[io + ix] - x_center)^2.0 + (y[jo + jy] - y_center)^2.0 <=
               (mountainwidth * width_factor)^2.0
                for iwm in 0:(spectral_modes - 1)
                    k_spectrum[iwm + 1, ix, jy] =
                        mountainwavenumber * cos(pi / spectral_modes * iwm)
                    l_spectrum[iwm + 1, ix, jy] =
                        mountainwavenumber * sin(pi / spectral_modes * iwm)
                    topography_spectrum[iwm + 1, ix, jy] =
                        0.5 *
                        mountainheight *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor * sqrt(
                                    (x[io + ix] - x_center)^2.0 +
                                    (y[jo + jy] - y_center)^2.0,
                                ),
                            )
                        ) / spectral_modes / (height_factor + 1.0)
                end
                topography_surface[ix, jy] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor * sqrt(
                                (x[io + ix] - x_center)^2.0 +
                                (y[jo + jy] - y_center)^2.0,
                            ),
                        )
                    ) *
                    height_factor / (height_factor + 1.0)
            end
        end
    end

    set_zonal_boundaries_of_field!(topography_surface, namelists, domain)
    set_meridional_boundaries_of_field!(topography_surface, namelists, domain)

    return (topography_surface, topography_spectrum, k_spectrum, l_spectrum)
end

"""
```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::MountainWave,
)
```

Compute topography for standard mountain wave test cases.

# Arguments

  - `namelists::Namelists`: Configuration parameters including domain and grid settings
  - `constants::Constants`: Physical constants including reference length scale
  - `domain::Domain`: Computational domain specification
  - `x::AbstractVector{<:AbstractFloat}`: X-coordinate grid points
  - `y::AbstractVector{<:AbstractFloat}`: Y-coordinate grid points
  - `testcase::MountainWave`: Standard mountain wave test case specification

# Returns

  - `Tuple`: (topography_surface, topography_spectrum, k_spectrum, l_spectrum)

      + `topography_surface`: 2D array of surface height values
      + `topography_spectrum`: Empty 3D array (not used for standard cases)
      + `k_spectrum`: Empty 3D array (not used for standard cases)
      + `l_spectrum`: Empty 3D array (not used for standard cases)

# Mountain Cases (Complete List)

## 2D Cases

  - **Case 1**: 2D cosine mountains

      + Formula: `h = 0.5 * H * (1 + cos(k * (x - x_c)))`
      + Simple periodic mountain ridge in x-direction

  - **Case 3**: 2D isolated mountain

      + Formula: `h = H / (1 + (x - x_c)²/W²)`
      + Classic bell-shaped isolated peak
  - **Case 5**: 2D cosine envelope with even background

      + Within envelope: Modulated cosine wave
      + Outside envelope: Constant background height
      + Smooth transition between regions
  - **Case 7**: 2D Gaussian envelope with even background

      + Formula: `h = 0.5 * H * (1 + exp(-(x-x_c)²/W²) * cos(k*(x-x_c)))`
      + Gaussian-modulated oscillatory structure
  - **Case 9**: 2D cosine envelope with cosine background

      + Within envelope: `h = 0.25 * H * envelope * (1 + cos(k*(x-x_c)))`
      + Localized cosine mountain with smooth edges
  - **Case 11**: 2D Gaussian envelope with Gaussian background

      + Formula: `h = 0.5 * H * exp(-(x-x_c)²/W²) * (1 + cos(k*(x-x_c)))`
      + Fully Gaussian-modulated structure

## 3D Cases

  - **Case 2**: 3D cosine mountains

      + Formula: `h = 0.5 * H * (1 + cos(k * √((x-x_c)² + (y-y_c)²)))`
      + Radially symmetric cosine structure

  - **Case 4**: 3D isolated mountain

      + Formula: `h = H / (1 + ((x-x_c)² + (y-y_c)²)/W²)`
      + Axisymmetric bell-shaped mountain
  - **Case 6**: 3D cosine envelope with even background

      + Within circular envelope: Modulated 3D cosine structure
      + Outside envelope: Constant background height
  - **Case 8**: 3D Gaussian envelope with even background

      + Formula: `h = 0.5 * H * (1 + exp(-r²/W²) * cos(k*r))`
      + Where `r = √((x-x_c)² + (y-y_c)²)`
  - **Case 10**: 3D cosine envelope with cosine background

      + Within envelope: Double cosine modulation in 3D
      + Localized 3D structure with smooth transitions
  - **Case 12**: 3D Gaussian envelope with Gaussian background

      + Formula: `h = 0.5 * H * exp(-r²/W²) * (1 + cos(k*r))`
      + Fully 3D Gaussian-modulated structure
  - **Case 13**: 3D WKB topography with spectral modes

      + Combines base topography with multiple spectral components
      + Each mode has specific wavenumber orientation
      + Used for complex wave interaction studies

# Parameters

  - `H`: Mountain height (mountainheight_dim/lref)
  - `W`: Mountain width (mountainwidth_dim/lref)
  - `k`: Mountain wavenumber (π/W)
  - `(x_c, y_c)`: Mountain center coordinates
  - `width_factor`: Envelope width scaling factor
  - `height_factor`: Height distribution factor (Case 13)
  - `spectral_modes`: Number of spectral components (Case 13)

# Notes

  - All coordinates are normalized by the reference length scale (lref)
  - Boundary conditions are applied to ensure proper domain periodicity
"""
function compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::AbstractTestCase,
)
    (; lx_dim, ly_dim) = namelists.domain
    (; testcase) = namelists.setting
    (;
        mountainheight_dim,
        mountainwidth_dim,
        mountain_case,
        height_factor,
        width_factor,
        spectral_modes,
    ) = namelists.grid
    (; nxx, nyy, io, jo, i0, i1, j0, j1) = domain
    (; lref) = constants

    mountainheight = mountainheight_dim / lref
    mountainwidth = mountainwidth_dim / lref
    mountainwavenumber = pi / mountainwidth

    x_center = 0.5 * (lx_dim[2] + lx_dim[1]) / lref
    y_center = 0.5 * (ly_dim[2] + ly_dim[1]) / lref

    topography_surface = zeros(nxx, nyy)
    topography_spectrum = zeros(0, 0, 0)
    k_spectrum = zeros(0, 0, 0)
    l_spectrum = zeros(0, 0, 0)

    for jy in j0:j1, ix in i0:i1
        if mountain_case == 1
            # 2D cosine mountains
            topography_surface[ix, jy] =
                0.5 *
                mountainheight *
                (1.0 + cos(mountainwavenumber * (x[io + ix] - x_center)))

        elseif mountain_case == 2
            # 3D cosine mountains
            topography_surface[ix, jy] =
                0.5 *
                mountainheight *
                (
                    1.0 + cos(
                        mountainwavenumber * sqrt(
                            (x[io + ix] - x_center)^2.0 +
                            (y[jo + jy] - y_center)^2.0,
                        ),
                    )
                )

        elseif mountain_case == 3
            # 2D isolated mountain
            topography_surface[ix, jy] =
                mountainheight /
                (1.0 + (x[io + ix] - x_center)^2.0 / mountainwidth^2.0)

        elseif mountain_case == 4
            # 3D isolated mountain
            topography_surface[ix, jy] =
                mountainheight / (
                    1.0 +
                    (
                        (x[io + ix] - x_center)^2.0 +
                        (y[jo + jy] - y_center)^2.0
                    ) / mountainwidth^2.0
                )

        elseif mountain_case == 5
            # 2D cosine envelope and even background
            if abs(x[io + ix] - x_center) <= mountainwidth * width_factor
                topography_surface[ix, jy] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 +
                        0.5 *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor *
                                (x[io + ix] - x_center),
                            )
                        ) *
                        cos(mountainwavenumber * (x[io + ix] - x_center))
                    )
            else
                topography_surface[ix, jy] = 0.5 * mountainheight
            end

        elseif mountain_case == 6
            # 3D cosine envelope and even background
            if (x[io + ix] - x_center)^2.0 + (y[jo + jy] - y_center)^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[ix, jy] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 +
                        0.5 *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor * sqrt(
                                    (x[io + ix] - x_center)^2.0 +
                                    (y[jo + jy] - y_center)^2.0,
                                ),
                            )
                        ) *
                        cos(
                            mountainwavenumber * sqrt(
                                (x[io + ix] - x_center)^2.0 +
                                (y[jo + jy] - y_center)^2.0,
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
                        -(x[io + ix] - x_center)^2.0 /
                        (mountainwidth * width_factor)^2.0,
                    ) * cos(mountainwavenumber * (x[io + ix] - x_center))
                )

        elseif mountain_case == 8
            # 3D Gaussian envelope and even background
            topography_surface[ix, jy] =
                0.5 *
                mountainheight *
                (
                    1.0 +
                    exp(
                        -(
                            (x[io + ix] - x_center)^2.0 +
                            (y[jo + jy] - y_center)^2.0
                        ) / (mountainwidth * width_factor)^2.0,
                    ) * cos(
                        mountainwavenumber * sqrt(
                            (x[io + ix] - x_center)^2.0 +
                            (y[jo + jy] - y_center)^2.0,
                        ),
                    )
                )

        elseif mountain_case == 9
            # 2D cosine envelope and cosine background
            if abs(x[io + ix] - x_center) <= mountainwidth * width_factor
                topography_surface[ix, jy] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            (x[io + ix] - x_center),
                        )
                    ) *
                    (1.0 + cos(mountainwavenumber * (x[io + ix] - x_center)))
            end

        elseif mountain_case == 10
            # 3D cosine envelope and cosine background
            if (x[io + ix] - x_center)^2.0 + (y[jo + jy] - y_center)^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[ix, jy] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor * sqrt(
                                (x[io + ix] - x_center)^2.0 +
                                (y[jo + jy] - y_center)^2.0,
                            ),
                        )
                    ) *
                    (
                        1.0 + cos(
                            mountainwavenumber * sqrt(
                                (x[io + ix] - x_center)^2.0 +
                                (y[jo + jy] - y_center)^2.0,
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
                    -(x[io + ix] - x_center)^2.0 /
                    (mountainwidth * width_factor)^2.0,
                ) *
                (1.0 + cos(mountainwavenumber * (x[io + ix] - x_center)))

        elseif mountain_case == 12
            # 3D Gaussian envelope and Gaussian background
            topography_surface[ix, jy] =
                0.5 *
                mountainheight *
                exp(
                    -(
                        (x[io + ix] - x_center)^2.0 +
                        (y[jo + jy] - y_center)^2.0
                    ) / (mountainwidth * width_factor)^2.0,
                ) *
                (
                    1.0 + cos(
                        mountainwavenumber * sqrt(
                            (x[io + ix] - x_center)^2.0 +
                            (y[jo + jy] - y_center)^2.0,
                        ),
                    )
                )

        elseif mountain_case == 13
            # 3D WKB topography
            if (x[io + ix] - x_center)^2.0 + (y[jo + jy] - y_center)^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[ix, jy] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor * sqrt(
                                (x[io + ix] - x_center)^2.0 +
                                (y[jo + jy] - y_center)^2.0,
                            ),
                        )
                    ) *
                    height_factor / (height_factor + 1.0)
                for iwm in 0:(spectral_modes - 1)
                    kk = mountainwavenumber * cos(pi / spectral_modes * iwm)
                    ll = mountainwavenumber * sin(pi / spectral_modes * iwm)
                    topography_surface[ix, jy] =
                        topography_surface[ix, jy] +
                        0.5 *
                        mountainheight *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor * sqrt(
                                    (x[io + ix] - x_center)^2.0 +
                                    (y[jo + jy] - y_center)^2.0,
                                ),
                            )
                        ) *
                        cos(
                            kk * (x[io + ix] - x_center) +
                            ll * (y[jo + jy] - y_center),
                        ) / spectral_modes / (height_factor + 1.0)
                end
            end
        else
            error("Error in compute_topography: Unknown mountain case!")
        end
    end

    set_zonal_boundaries_of_field!(topography_surface, namelists, domain)
    set_meridional_boundaries_of_field!(topography_surface, namelists, domain)

    return (topography_surface, topography_spectrum, k_spectrum, l_spectrum)
end

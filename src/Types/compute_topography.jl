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
    (; nxx, nyy, io, jo) = domain
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

    for jy in 1:nyy, ix in 1:nxx
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

        else
        end
    end

    return (topography_surface, topography_spectrum, k_spectrum, l_spectrum)
end

function compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::MountainWave,
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
    (; nxx, nyy, io, jo) =
        domain
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

    for jy in 1:nyy, ix in 1:nxx
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

    return (topography_surface, topography_spectrum, k_spectrum, l_spectrum)
end
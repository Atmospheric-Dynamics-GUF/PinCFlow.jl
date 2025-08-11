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

Compute and return the topography for the WKB-mountain-wave test case.

The supported topography shapes are as follows, listed according to the value of `namelists.grid.mountain_case`.

 1. 2D cosine mountains:

    ```math
        h_\\mathrm{b} = \\frac{h_0}{2}, \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0, \\quad h_\\mathrm{w} = \\frac{h_0}{2}
    ```

 1. ``-``

 1. ``-``

 1. ``-``

 1. 2D cosine envelope and even background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} & = \\frac{h_0}{2}, \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\begin{cases}
            \\frac{h_0}{4} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\left(x - x_0\\right)\\right]\\right\\} & \\mathrm{if} \\quad \\left|x - x_0\\right| \\leq r_l l_0,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

 1. ``-``

 1. 2D Gaussian envelope and even background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} & = \\frac{h_0}{2}, \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\frac{h_0}{2} \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2\\right]
    \\end{align*}
    ```

 1. ``-``

 1. 2D cosine envelope and cosine background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x\\right) & = h_\\mathrm{w} \\left(x\\right), \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\begin{cases}
            \\frac{h_0}{4} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\left(x - x_0\\right)\\right]\\right\\} & \\mathrm{if} \\quad \\left|x - x_0\\right| \\leq r_l l_0,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

 1. ``-``

 1. 2D Gaussian envelope and Gaussian background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x\\right) & = h_\\mathrm{w} \\left(x\\right), \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\frac{h_0}{2} \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2\\right]
    \\end{align*}
    ```

 1. ``-``

 1. 3D WKB topography:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x, y\\right) & = r_h n_h h_\\mathrm{w} \\left(x, y\\right), \\quad k_{h, \\alpha} = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_{h, \\alpha} = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right),\\\\
        h_\\mathrm{w} \\left(x, y\\right) & = \\begin{cases}
            \\frac{h_0}{2 n_h \\left(r_h + 1\\right)} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\} & \\mathrm{if} \\quad \\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2 \\leq r_l^2 l_0^2,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

Therein, ``h_0``, ``l_0``, ``r_h``, ``r_l`` and ``n_h`` are given by the properties `mountainheight_dim`, `mountainwidth_dim`, `height_factor`, `width_factor` and `spectral_modes` of `namelists.grid`, respectively, whereas ``\\left(x_0, y_0\\right)`` is the horizontal center of the domain.

```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::AbstractTestCase,
)
```

Compute and return the topography for non-WKB-mountain-wave test cases.

The supported topography shapes are as follows, listed according to the value of `namelists.grid.mountain_case`.

 1. 2D cosine mountains:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\left(x - x_0\\right)\\right]\\right\\}
    ```

 1. 3D cosine mountains:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\}
    ```

 1. 2D isolated mountain:

    ```math
    h \\left(x\\right) = \\frac{h_0}{1 + \\left(x - x_0\\right)^2 / l_0^2}
    ```

 1. 3D isolated mountain:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{1 + \\left[\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2\\right] / l_0^2}
    ```

 1. 2D cosine envelope and even background:

    ```math
    h \\left(x\\right) = \\begin{cases}
        \\frac{h_0}{2} \\left\\{1 + \\frac{1}{2} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\left(x - x_0\\right)\\right)\\right] \\cos \\left[\\frac{\\pi}{l_0} \\left(x - x_0\\right)\\right]\\right\\} & \\mathrm{if} \\quad \\left|x - x_0\\right| \\leq r_l l_0,\\\\
        \\frac{h_0}{2} & \\mathrm{else}
    \\end{cases}
    ```

 1. 3D cosine envelope and even background:

    ```math
    h \\left(x, y\\right) = \\begin{cases}
        \\frac{h_0}{2} \\left\\{1 + \\frac{1}{2} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right)\\right] \\cos \\left[\\frac{\\pi}{l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\} & \\mathrm{if} \\quad \\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2 \\leq r_l^2 l_0^2,\\\\
        \\frac{h_0}{2} & \\mathrm{else}
    \\end{cases}
    ```

 1. 2D Gaussian envelope and even background:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\left\\{1 + \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2\\right] \\cos \\left[\\frac{\\pi}{l_0} \\left(x - x_0\\right)\\right]\\right\\}
    ```

 1. 3D Gaussian envelope and even background:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\left\\{1 + \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2 - \\left(\\frac{y - y_0}{r_l l_0}\\right)^2\\right] \\cos \\left[\\frac{\\pi}{l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\}
    ```

 1. 2D cosine envelope and cosine background:

    ```math
    h \\left(x\\right) = \\begin{cases}
        \\frac{h_0}{4} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\left(x - x_0\\right)\\right]\\right\\} \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\left(x - x_0\\right)\\right]\\right\\} & \\mathrm{if} \\quad \\left|x - x_0\\right| \\leq r_l l_0,\\\\
        0 & \\mathrm{else}
    \\end{cases}
    ```

 1. 3D cosine envelope and cosine background:

    ```math
    h \\left(x, y\\right) = \\begin{cases}
        \\frac{h_0}{4} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\} \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\} & \\mathrm{if} \\quad \\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2 \\leq r_l^2 l_0^2,\\\\
        0 & \\mathrm{else}
    \\end{cases}
    ```

 1. 2D Gaussian envelope and Gaussian background:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2\\right] \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\left(x - x_0\\right)\\right]\\right\\}
    ```

 1. 3D Gaussian envelope and Gaussian background:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\exp \\left[- \\left(\\frac{x - x_0}{r_l l_0}\\right)^2 - \\left(\\frac{y - y_0}{r_l l_0}\\right)^2\\right] \\left\\{1 + \\cos \\left[\\frac{\\pi}{l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\}
    ```

 1. 3D WKB topography:

    ```math
    \\begin{align*}
        h \\left(x, y\\right) & = \\begin{cases}
            \\frac{h_0}{2 \\left(r_h + 1\\right)} \\left\\{1 + \\cos \\left[\\frac{\\pi}{r_l l_0} \\sqrt{\\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2}\\right]\\right\\} \\left\\{r_h + n_h^{- 1} \\sum\\limits_{\\alpha = 0}^{n_h - 1} \\cos \\left[k_\\alpha \\left(x - x_0\\right) + l_\\alpha \\left(y - y_0\\right)\\right]\\right\\} & \\mathrm{if} \\quad \\left(x - x_0\\right)^2 + \\left(y - y_0\\right)^2 \\leq r_l^2 l_0^2,\\\\
            0 & \\mathrm{else},
        \\end{cases}\\\\
        k_\\alpha & = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_\\alpha = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right)
    \\end{align*}
    ```

Therein, ``h_0``, ``l_0``, ``r_h``, ``r_l`` and ``n_h`` are given by the properties `mountainheight_dim`, `mountainwidth_dim`, `height_factor`, `width_factor` and `spectral_modes` of `namelists.grid`, respectively, whereas ``\\left(x_0, y_0\\right)`` is the horizontal center of the domain. The arrays representing the unresolved spectrum are set to have the size `(0, 0, 0)`.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `x`: ``\\widehat{x}``-coordinate grid points.
  - `y`: ``\\widehat{y}``-coordinate grid points.
  - `testcase`: Test case on which the current simulation is based.

# Returns

  - `::AbstractMatrix{<:AbstractFloat}`: Resolved orography.
  - `::AbstractArray{<:AbstractFloat, 3}`: Spectrum of the unresolved orography.
  - `::AbstractArray{<:AbstractFloat, 3}`: Zonal wavenumbers of the orographic spectrum.
  - `::AbstractArray{<:AbstractFloat, 3}`: Meridional wavenumbers of the orographic spectrum.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)
  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_topography end

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

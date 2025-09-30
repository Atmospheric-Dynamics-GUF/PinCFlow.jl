"""
```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::AbstractWKBTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
```

Compute and return the topography for WKB test cases.

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
            \\frac{h_0}{4} \\left[1 + \\cos \\left(\\frac{\\pi x}{r_l l_0}\\right)\\right] & \\mathrm{if} \\quad \\left|x\\right| \\leq r_l l_0,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

 1. ``-``

 1. 2D Gaussian envelope and even background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} & = \\frac{h_0}{2}, \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\frac{h_0}{2} \\exp \\left(- \\frac{x^2}{r_l^2 l_0^2}\\right)
    \\end{align*}
    ```

 1. ``-``

 1. 2D cosine envelope and cosine background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x\\right) & = h_\\mathrm{w} \\left(x\\right), \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\begin{cases}
            \\frac{h_0}{4} \\left[1 + \\cos \\left(\\frac{\\pi x}{r_l l_0}\\right)\\right] & \\mathrm{if} \\quad \\left|x\\right| \\leq r_l l_0,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

 1. ``-``

 1. 2D Gaussian envelope and Gaussian background:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x\\right) & = h_\\mathrm{w} \\left(x\\right), \\quad k_h = \\frac{\\pi}{l_0}, \\quad l_h = 0,\\\\
        h_\\mathrm{w} \\left(x\\right) & = \\frac{h_0}{2} \\exp \\left(- \\frac{x^2}{r_l^2 l_0^2}\\right)
    \\end{align*}
    ```

 1. ``-``

 1. 3D WKB topography:

    ```math
    \\begin{align*}
        h_\\mathrm{b} \\left(x, y\\right) & = r_h n_h h_\\mathrm{w} \\left(x, y\\right), \\quad k_{h, \\alpha} = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_{h, \\alpha} = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right),\\\\
        h_\\mathrm{w} \\left(x, y\\right) & = \\begin{cases}
            \\frac{h_0}{2 n_h \\left(r_h + 1\\right)} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
            0 & \\mathrm{else}
        \\end{cases}
    \\end{align*}
    ```

Therein, ``h_0``, ``l_0``, ``r_h``, ``r_l`` and ``n_h`` are given by the properties `mountain_height`, `mountain_half_width`, `height_factor`, `width_factor` and `spectral_modes` of `namelists.grid`, respectively.

The arrays in the returned tuple represent (in order) the resolved topography, the amplitudes of the unresolved topography, the corresponding zonal wavenumbers and the corresponding meridional wavenumbers.

```julia
compute_topography(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    x::AbstractVector{<:AbstractFloat},
    y::AbstractVector{<:AbstractFloat},
    testcase::AbstractTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
```

Compute and return the topography for non-WKB test cases.

The supported topography shapes are as follows, listed according to the value of `namelists.grid.mountain_case`.

 1. 2D cosine mountains:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\left[1 + \\cos \\left(\\frac{\\pi x}{l_0}\\right)\\right]
    ```

 1. 3D cosine mountains:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\left[1 + \\cos \\left(\\frac{\\pi}{l_0} \\sqrt{x^2 + y^2}\\right)\\right]
    ```

 1. 2D isolated mountain:

    ```math
    h \\left(x\\right) = \\frac{h_0}{1 + x^2 / l_0^2}
    ```

 1. 3D isolated mountain:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{1 + \\left(x^2 + y^2\\right) / l_0^2}
    ```

 1. 2D cosine envelope and even background:

    ```math
    h \\left(x\\right) = \\begin{cases}
        \\frac{h_0}{2} \\left\\{1 + \\frac{1}{2} \\left[1 + \\cos \\left(\\frac{\\pi x}{r_l l_0}\\right)\\right] \\cos \\left(\\frac{\\pi x}{l_0}\\right)\\right\\} & \\mathrm{if} \\quad \\left|x\\right| \\leq r_l l_0,\\\\
        \\frac{h_0}{2} & \\mathrm{else}
    \\end{cases}
    ```

 1. 3D cosine envelope and even background:

    ```math
    h \\left(x, y\\right) = \\begin{cases}
        \\frac{h_0}{2} \\left\\{1 + \\frac{1}{2} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] \\cos \\left(\\frac{\\pi}{l_0} \\sqrt{x^2 + y^2}\\right)\\right\\} & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
        \\frac{h_0}{2} & \\mathrm{else}
    \\end{cases}
    ```

 1. 2D Gaussian envelope and even background:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\left[1 + \\exp \\left(- \\frac{x^2}{r_l^2 l_0^2}\\right) \\cos \\left(\\frac{\\pi x}{l_0}\\right)\\right]
    ```

 1. 3D Gaussian envelope and even background:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\left[1 + \\exp \\left(- \\frac{x^2 + y^2}{r_l^2 l_0^2}\\right) \\cos \\left(\\frac{\\pi}{l_0} \\sqrt{x^2 + y^2}\\right)\\right]
    ```

 1. 2D cosine envelope and cosine background:

    ```math
    h \\left(x\\right) = \\begin{cases}
        \\frac{h_0}{4} \\left[1 + \\cos \\left(\\frac{\\pi x}{r_l l_0}\\right)\\right] \\left[1 + \\cos \\left(\\frac{\\pi x}{l_0}\\right)\\right] & \\mathrm{if} \\quad \\left|x\\right| \\leq r_l l_0,\\\\
        0 & \\mathrm{else}
    \\end{cases}
    ```

 1. 3D cosine envelope and cosine background:

    ```math
    h \\left(x, y\\right) = \\begin{cases}
        \\frac{h_0}{4} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] \\left[1 + \\cos \\left(\\frac{\\pi}{l_0} \\sqrt{x^2 + y^2}\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
        0 & \\mathrm{else}
    \\end{cases}
    ```

 1. 2D Gaussian envelope and Gaussian background:

    ```math
    h \\left(x\\right) = \\frac{h_0}{2} \\exp \\left(- \\frac{x^2}{r_l^2 l_0^2}\\right) \\left[1 + \\cos \\left(\\frac{\\pi x}{l_0}\\right)\\right]
    ```

 1. 3D Gaussian envelope and Gaussian background:

    ```math
    h \\left(x, y\\right) = \\frac{h_0}{2} \\exp \\left(- \\frac{x^2 + y^2}{r_l^2 l_0^2}\\right) \\left[1 + \\cos \\left(\\frac{\\pi}{l_0} \\sqrt{x^2 + y^2}\\right)\\right]
    ```

 1. 3D WKB topography:

    ```math
    {\\small\\begin{align*}
        h \\left(x, y\\right) & = \\begin{cases}
            \\frac{h_0}{2 \\left(r_h + 1\\right)} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] \\left[r_h + n_h^{- 1} \\sum\\limits_{\\alpha = 0}^{n_h - 1} \\cos \\left(k_\\alpha x + l_\\alpha y\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
            0 & \\mathrm{else},
        \\end{cases}\\\\
        k_\\alpha & = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_\\alpha = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right)
    \\end{align*}}
    ```

Therein, ``h_0``, ``l_0``, ``r_h``, ``r_l`` and ``n_h`` are given by the properties `mountain_height`, `mountain_half_width`, `height_factor`, `width_factor` and `spectral_modes` of `namelists.grid`, respectively. The arrays representing the unresolved spectrum are set to have the size `(0, 0, 0)`.

The topography is represented by the first array in the returned tuple.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `x`: ``\\widehat{x}``-coordinate grid points.

  - `y`: ``\\widehat{y}``-coordinate grid points.

  - `testcase`: Test case on which the current simulation is based.

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
    testcase::AbstractWKBTestCase,
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
    (; testcase) = namelists.setting
    (; nwm) = namelists.wkb
    (;
        mountain_height,
        mountain_half_width,
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

    mountainheight = mountain_height / lref
    mountainwidth = mountain_half_width / lref
    mountainwavenumber = pi / mountainwidth

    topography_surface = zeros(nxx, nyy)
    topography_spectrum = zeros(nwm, nxx, nyy)
    k_spectrum = zeros(nwm, nxx, nyy)
    l_spectrum = zeros(nwm, nxx, nyy)

    @ivy for j in j0:j1, i in i0:i1
        if mountain_case == 1
            # 2D cosine mountains
            topography_surface[i, j] = 0.5 * mountainheight
            topography_spectrum[1, i, j] = 0.5 * mountainheight
            k_spectrum[1, i, j] = mountainwavenumber

        elseif mountain_case == 5
            # 2D cosine envelope and even background
            if abs(x[io + i]) <= mountainwidth * width_factor
                k_spectrum[1, i, j] = mountainwavenumber
                topography_spectrum[1, i, j] =
                    0.25 *
                    mountainheight *
                    (1.0 + cos(mountainwavenumber / width_factor * x[io + i]))
            end
            topography_surface[i, j] = 0.5 * mountainheight

        elseif mountain_case == 7
            # 2D Gaussian envelope and even background
            k_spectrum[1, i, j] = mountainwavenumber
            topography_spectrum[1, i, j] =
                0.5 *
                mountainheight *
                exp(-x[io + i]^2.0 / (mountainwidth * width_factor)^2.0)
            topography_surface[i, j] = 0.5 * mountainheight

        elseif mountain_case == 9
            # 2D cosine envelope and cosine background
            if abs(x[io + i]) <= mountainwidth * width_factor
                k_spectrum[1, i, j] = mountainwavenumber
                topography_spectrum[1, i, j] =
                    0.25 *
                    mountainheight *
                    (1.0 + cos(mountainwavenumber / width_factor * x[io + i]))
                topography_surface[i, j] =
                    0.25 *
                    mountainheight *
                    (1.0 + cos(mountainwavenumber / width_factor * x[io + i]))
            end

        elseif mountain_case == 11
            # 2D Gaussian envelope and Gaussian background
            k_spectrum[1, i, j] = mountainwavenumber
            topography_spectrum[1, i, j] =
                0.5 *
                mountainheight *
                exp(-x[io + i]^2.0 / (mountainwidth * width_factor)^2.0)
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                exp(-x[io + i]^2.0 / (mountainwidth * width_factor)^2.0)

        elseif mountain_case == 13
            # 3D WKB topography
            if x[io + i]^2.0 + y[jo + j]^2.0 <=
               (mountainwidth * width_factor)^2.0
                for alpha in 0:(spectral_modes - 1)
                    k_spectrum[alpha + 1, i, j] =
                        mountainwavenumber * cos(pi / spectral_modes * alpha)
                    l_spectrum[alpha + 1, i, j] =
                        mountainwavenumber * sin(pi / spectral_modes * alpha)
                    topography_spectrum[alpha + 1, i, j] =
                        0.5 *
                        mountainheight *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor *
                                sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                            )
                        ) / spectral_modes / (height_factor + 1.0)
                end
                topography_surface[i, j] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
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
)::Tuple{
    <:AbstractMatrix{<:AbstractFloat},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
    <:AbstractArray{<:AbstractFloat, 3},
}
    (; testcase) = namelists.setting
    (;
        mountain_height,
        mountain_half_width,
        mountain_case,
        height_factor,
        width_factor,
        spectral_modes,
    ) = namelists.grid
    (; nxx, nyy, io, jo, i0, i1, j0, j1) = domain
    (; lref) = constants

    mountainheight = mountain_height / lref
    mountainwidth = mountain_half_width / lref
    mountainwavenumber = pi / mountainwidth

    topography_surface = zeros(nxx, nyy)
    topography_spectrum = zeros(0, 0, 0)
    k_spectrum = zeros(0, 0, 0)
    l_spectrum = zeros(0, 0, 0)

    @ivy for j in j0:j1, i in i0:i1
        if mountain_case == 1
            # 2D cosine mountains
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                (1.0 + cos(mountainwavenumber * x[io + i]))

        elseif mountain_case == 2
            # 3D cosine mountains
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                (
                    1.0 + cos(
                        mountainwavenumber *
                        sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                    )
                )

        elseif mountain_case == 3
            # 2D isolated mountain
            topography_surface[i, j] =
                mountainheight / (1.0 + x[io + i]^2.0 / mountainwidth^2.0)

        elseif mountain_case == 4
            # 3D isolated mountain
            topography_surface[i, j] =
                mountainheight /
                (1.0 + (x[io + i]^2.0 + y[jo + j]^2.0) / mountainwidth^2.0)

        elseif mountain_case == 5
            # 2D cosine envelope and even background
            if abs(x[io + i]) <= mountainwidth * width_factor
                topography_surface[i, j] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 +
                        0.5 *
                        (
                            1.0 +
                            cos(mountainwavenumber / width_factor * x[io + i])
                        ) *
                        cos(mountainwavenumber * x[io + i])
                    )
            else
                topography_surface[i, j] = 0.5 * mountainheight
            end

        elseif mountain_case == 6
            # 3D cosine envelope and even background
            if x[io + i]^2.0 + y[jo + j]^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[i, j] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 +
                        0.5 *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor *
                                sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                            )
                        ) *
                        cos(
                            mountainwavenumber *
                            sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                        )
                    )
            else
                topography_surface[i, j] = 0.5 * mountainheight
            end

        elseif mountain_case == 7
            # 2D Gaussian envelope and even background
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                (
                    1.0 +
                    exp(-x[io + i]^2.0 / (mountainwidth * width_factor)^2.0) *
                    cos(mountainwavenumber * x[io + i])
                )

        elseif mountain_case == 8
            # 3D Gaussian envelope and even background
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                (
                    1.0 +
                    exp(
                        -(x[io + i]^2.0 + y[jo + j]^2.0) /
                        (mountainwidth * width_factor)^2.0,
                    ) * cos(
                        mountainwavenumber *
                        sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                    )
                )

        elseif mountain_case == 9
            # 2D cosine envelope and cosine background
            if abs(x[io + i]) <= mountainwidth * width_factor
                topography_surface[i, j] =
                    0.25 *
                    mountainheight *
                    (1.0 + cos(mountainwavenumber / width_factor * x[io + i])) *
                    (1.0 + cos(mountainwavenumber * x[io + i]))
            end

        elseif mountain_case == 10
            # 3D cosine envelope and cosine background
            if x[io + i]^2.0 + y[jo + j]^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[i, j] =
                    0.25 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                        )
                    ) *
                    (
                        1.0 + cos(
                            mountainwavenumber *
                            sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                        )
                    )
            end

        elseif mountain_case == 11
            # 2D Gaussian envelope and Gaussian background
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                exp(-x[io + i]^2.0 / (mountainwidth * width_factor)^2.0) *
                (1.0 + cos(mountainwavenumber * x[io + i]))

        elseif mountain_case == 12
            # 3D Gaussian envelope and Gaussian background
            topography_surface[i, j] =
                0.5 *
                mountainheight *
                exp(
                    -(x[io + i]^2.0 + y[jo + j]^2.0) /
                    (mountainwidth * width_factor)^2.0,
                ) *
                (
                    1.0 + cos(
                        mountainwavenumber *
                        sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                    )
                )

        elseif mountain_case == 13
            # 3D WKB topography
            if x[io + i]^2.0 + y[jo + j]^2.0 <=
               (mountainwidth * width_factor)^2.0
                topography_surface[i, j] =
                    0.5 *
                    mountainheight *
                    (
                        1.0 + cos(
                            mountainwavenumber / width_factor *
                            sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                        )
                    ) *
                    height_factor / (height_factor + 1.0)
                for alpha in 0:(spectral_modes - 1)
                    kk = mountainwavenumber * cos(pi / spectral_modes * alpha)
                    ll = mountainwavenumber * sin(pi / spectral_modes * alpha)
                    topography_surface[i, j] =
                        topography_surface[i, j] +
                        0.5 *
                        mountainheight *
                        (
                            1.0 + cos(
                                mountainwavenumber / width_factor *
                                sqrt(x[io + i]^2.0 + y[jo + j]^2.0),
                            )
                        ) *
                        cos(kk * x[io + i] + ll * y[jo + j]) / spectral_modes /
                        (height_factor + 1.0)
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

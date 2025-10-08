"""
```julia
mountain_range(
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
```

Return the resolved part of the configured mountain range at ``\\left(x, y\\right)`` by dispatching to a test-case-specific method.

```julia
mountain_range(
    test_case::AbstractTestCase,
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
```

Return the configured mountain-range at ``\\left(x, y\\right)`` for non-WKB test cases.

The mountain-range is given by

```math
{\\small\\begin{align*}
    h \\left(x, y\\right) & = \\begin{cases}
        \\frac{h_0}{2 \\left(r_h + 1\\right)} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] \\left[r_h + n_h^{- 1} \\sum\\limits_{\\alpha = 0}^{n_h - 1} \\cos \\left(k_\\alpha x + l_\\alpha y\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
        0 & \\mathrm{else},
    \\end{cases}\\\\
    k_\\alpha & = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_\\alpha = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right),
\\end{align*}}
```

where ``h_0``, ``l_0``, ``r_h``, ``r_l`` and ``n_h`` are given by the properties `mountain_height`, `mountain_half_width`, `height_factor`, `width_factor` and `spectral_modes` of `namelists.grid`, respectively.

```julia
mountain_range(
    test_case::AbstractWKBTestCase,
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
```

Return the resolved part of the configured mountain-range at ``\\left(x, y\\right)`` for WKB test cases.

The resolved part is given by

```math
\\begin{align*}
    h_\\mathrm{b} \\left(x, y\\right) & = \\begin{cases}
        \\frac{r_h h_0}{2 \\left(r_h + 1\\right)} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
        0 & \\mathrm{else}.
    \\end{cases}
\\end{align*}
```

```julia
mountain_range(
    namelists::Namelists,
    alpha::Integer,
    x::AbstractFloat,
    y::AbstractFloat,
)::NTuple{3, <:AbstractFloat}
```

Return the unresolved part of the configured mountain-range at ``\\left(x, y\\right)`` for WKB test cases (i.e. zonal wavenumber, meridional wavenumber and wave ampltitude).

The unresolved part is given by

```math
\\begin{align*}
    k_{h, \\alpha} & = \\frac{\\pi}{l_0} \\cos \\left(\\frac{\\pi \\alpha}{n_h}\\right), \\quad l_{h, \\alpha} = \\frac{\\pi}{l_0} \\sin \\left(\\frac{\\pi \\alpha}{n_h}\\right),\\\\
    h_{\\mathrm{w}, \\alpha} \\left(x, y\\right) & = \\begin{cases}
        \\frac{h_0}{2 n_h \\left(r_h + 1\\right)} \\left[1 + \\cos \\left(\\frac{\\pi}{r_l l_0} \\sqrt{x^2 + y^2}\\right)\\right] & \\mathrm{if} \\quad x^2 + y^2 \\leq r_l^2 l_0^2,\\\\
        0 & \\mathrm{else}.
    \\end{cases}
\\end{align*}
```

# Arguments

  - `test_case`: Test case on which the current simulation is based.

  - `namelists`: Namelists with all model parameters.

  - `alpha`: Wave-mode index.

  - `x`: Zonal position.

  - `y`: Meridional position.
"""
function mountain_range end

function mountain_range(
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
    (; test_case) = namelists.setting
    return mountain_range(test_case, namelists, x, y)
end

function mountain_range(
    test_case::AbstractTestCase,
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
    (;
        mountain_height,
        mountain_half_width,
        height_factor,
        width_factor,
        spectral_modes,
    ) = namelists.grid

    h0 = mountain_height
    l0 = mountain_half_width
    rh = height_factor
    rl = width_factor
    nh = spectral_modes

    if x^2 + y^2 <= (l0 * rl)^2
        hb =
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh / (rh + 1)
        for alpha in 0:(nh - 1)
            k = pi / l0 * cos(pi / nh * alpha)
            l = pi / l0 * sin(pi / nh * alpha)
            hb +=
                h0 / 2 *
                (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) *
                cos(k * x + l * y) / nh / (rh + 1.0)
        end
    else
        hb = 0.0
    end

    return hb
end

function mountain_range(
    test_case::AbstractWKBTestCase,
    namelists::Namelists,
    x::AbstractFloat,
    y::AbstractFloat,
)::AbstractFloat
    (; mountain_height, mountain_half_width, height_factor, width_factor) =
        namelists.grid

    h0 = mountain_height
    l0 = mountain_half_width
    rh = height_factor
    rl = width_factor

    if x^2 + y^2 <= (rl * l0)^2
        hb =
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) * rh / (rh + 1)
    else
        hb = 0.0
    end

    return hb
end

function mountain_range(
    namelists::Namelists,
    alpha::Integer,
    x::AbstractFloat,
    y::AbstractFloat,
)::NTuple{3, <:AbstractFloat}
    (;
        mountain_height,
        mountain_half_width,
        height_factor,
        width_factor,
        spectral_modes,
    ) = namelists.grid

    h0 = mountain_height
    l0 = mountain_half_width
    rh = height_factor
    rl = width_factor
    nh = spectral_modes

    if x^2 + y^2 <= (rl * l0)^2
        kh = pi / l0 * cos(pi / nh * (alpha - 1))
        lh = pi / l0 * sin(pi / nh * (alpha - 1))
        hw =
            h0 / 2 * (1 + cos(pi / (rl * l0) * sqrt(x^2 + y^2))) / nh / (rh + 1)
    else
        kh = lh = hw = 0.0
    end

    return (kh, lh, hw)
end

"""
```julia
interpolate(
    namelists::Namelists;
    philbd::AbstractFloat = NaN,
    philbu::AbstractFloat = NaN,
    philfd::AbstractFloat = NaN,
    philfu::AbstractFloat = NaN,
    phirbd::AbstractFloat = NaN,
    phirbu::AbstractFloat = NaN,
    phirfd::AbstractFloat = NaN,
    phirfu::AbstractFloat = NaN,
    zlbd::AbstractFloat = NaN,
    zlbu::AbstractFloat = NaN,
    zlfd::AbstractFloat = NaN,
    zlfu::AbstractFloat = NaN,
    zrbd::AbstractFloat = NaN,
    zrbu::AbstractFloat = NaN,
    zrfd::AbstractFloat = NaN,
    zrfu::AbstractFloat = NaN,
    zlc::AbstractFloat = NaN,
    yb::AbstractFloat = NaN,
    yf::AbstractFloat = NaN,
    ylc::AbstractFloat = NaN,
    xl::AbstractFloat = NaN,
    xr::AbstractFloat = NaN,
    xlc::AbstractFloat = NaN,
)
```

Perform trilinear interpolation to `(xlc, ylc, zlc)`, with values from eight surrounding grid points (two zonal positions, two meridional positions and eight vertical positions).

Out of the eight grid points, four each are assumed to be to the left, to the right, behind, in front of, below and above the location of interest. Due to the grid being terrain-following, this includes eight different vertical positions, but only two zonal and two meridional positions. This is handled by performing successive linear interpolations, where the vertical position is interpolated along with the field of interest.

The exact algorithm is as follows.

  1. Interpolation in ``x``:

```math
\\begin{align*}
\\psi_\\mathrm{BD} & = f_x \\psi_\\mathrm{LBD} + (1 - f_x) \\psi_\\mathrm{RBD},\\\\
\\psi_\\mathrm{BU} & = f_x \\psi_\\mathrm{LBU} + (1 - f_x) \\psi_\\mathrm{RBU},\\\\
\\psi_\\mathrm{FD} & = f_x \\psi_\\mathrm{LFD} + (1 - f_x) \\psi_\\mathrm{RFD},\\\\
\\psi_\\mathrm{FU} & = f_x \\psi_\\mathrm{LFU} + (1 - f_x) \\psi_\\mathrm{RFU}
\\end{align*}
```

  2. Interpolation in ``y``:

```math
\\begin{align*}
\\psi_\\mathrm{D} & = f_y \\psi_\\mathrm{BD} + (1 - f_y) \\psi_\\mathrm{FD},\\\\
\\psi_\\mathrm{U} & = f_y \\psi_\\mathrm{BU} + (1 - f_y) \\psi_\\mathrm{FU}
\\end{align*}
```

  3. Interpolation in ``z``:

```math
\\phi_\\mathrm{C} = f_z \\phi_\\mathrm{D} + (1 - f_z) \\phi_\\mathrm{U}
```

Therein, the acronyms ``\\mathrm{L}``, ``\\mathrm{R}``, ``\\mathrm{B}``, ``\\mathrm{F}``, ``\\mathrm{D}`` and ``\\mathrm{U}`` represent the grid points to the left, to the right, forward, backward, downward and upward of the location of interest (denoted by ``\\mathrm{C}``), respectively, ``\\psi = \\left(\\phi, z\\right)`` and

```math
\\f_\\alpha = \\begin{cases}
0 & \\mathrm{if} \\quad \\alpha_\\beta = \\alpha_\\gamma,\\\\
1 & \\mathrm{if} \\quad \\alpha_\\mathrm{C} < \\alpha_\\beta,\\\\
\\frac{\\alpha_\\gamma - \\alpha_\\mathrm{C}}{\\alpha_\\gamma - \\alpha_\\beta} & \\alpha_\\beta \\leq \\alpha_\\mathrm{C} \\leq \\alpha_\\gamma,\\\\
0 & \\alpha_\\gamma < \\alpha_\\mathrm{C}
\\end{cases}
```

where ``\\left(\\alpha, \\beta, \\gamma\\right) \\in \\left\\{\\left(x, \\mathrm{L}, \\mathrm{R}\\right), \\left(y, \\mathrm{B}, \\mathrm{F}\\right), \\left(z, \\mathrm{D}, \\mathrm{U}\\right)\\right}``.

Due to their large number, the positions and values are given as keyword arguments with the default value `NaN`, so that their order does not matter and missing arguments are easy to detect.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `philbd`: Value at the point to the left, behind and below.
  - `philbu`: Value at the point to the left, behind and above.
  - `philfd`: Value at the point to the left, in front and below.
  - `philfu`: Value at the point to the left, in front and above.
  - `phirbd`: Value at the point to the right, behind and below.
  - `phirbu`: Value at the point to the right, behind and above.
  - `phirfd`: Value at the point to the right, in front and below.
  - `phirfu`: Value at the point to the right, in front and above.
  - `zlbd`: Vertical coordinate of the point to the left, behind and below.
  - `zlbu`: Vertical coordinate of the point to the left, behind and above.
  - `zlfd`: Vertical coordinate of the point to the left, in front and below.
  - `zlfu`: Vertical coordinate of the point to the left, in front and above.
  - `zrbd`: Vertical coordinate of the point to the right, behind and below.
  - `zrbu`: Vertical coordinate of the point to the right, behind and above.
  - `zrfd`: Vertical coordinate of the point to the right, in front and below.
  - `zrfu`: Vertical coordinate of the point to the right, in front and above.
  - `zlc`: Vertical position of interest.
  - `yb`: Meridional coordinate of the points behind.
  - `yf`: Meridional coordinate of the points in front.
  - `ylc`: Meridional position of interest.
  - `xl`: Zonal coordinate of the points to the left.
  - `xr`: Zonal coordinate of the points to the right.
  - `xlc`: Zonal position of interest.

# Returns

  - `::AbstractFloat`: Interpolated field value at the location of interest.
"""
function interpolate(
    namelists::Namelists;
    philbd::AbstractFloat = NaN,
    philbu::AbstractFloat = NaN,
    philfd::AbstractFloat = NaN,
    philfu::AbstractFloat = NaN,
    phirbd::AbstractFloat = NaN,
    phirbu::AbstractFloat = NaN,
    phirfd::AbstractFloat = NaN,
    phirfu::AbstractFloat = NaN,
    zlbd::AbstractFloat = NaN,
    zlbu::AbstractFloat = NaN,
    zlfd::AbstractFloat = NaN,
    zlfu::AbstractFloat = NaN,
    zrbd::AbstractFloat = NaN,
    zrbu::AbstractFloat = NaN,
    zrfd::AbstractFloat = NaN,
    zrfu::AbstractFloat = NaN,
    zlc::AbstractFloat = NaN,
    yb::AbstractFloat = NaN,
    yf::AbstractFloat = NaN,
    ylc::AbstractFloat = NaN,
    xl::AbstractFloat = NaN,
    xr::AbstractFloat = NaN,
    xlc::AbstractFloat = NaN,
)
    (; sizex, sizey) = namelists.domain

    # Interpolate in x.
    if sizex == 1
        phibd = philbd
        phibu = philbu

        phifd = philfd
        phifu = philfu

        zbd = zlbd
        zbu = zlbu

        zfd = zlfd
        zfu = zlfu
    else
        if xr < xl
            error("Error in interpolate: xr = ", xr, " < xl = ", xl)
        elseif xr == xl
            factor = 0.0
        elseif xlc > xr
            factor = 0.0
        elseif xlc > xl
            factor = (xr - xlc) / (xr - xl)
        else
            factor = 1.0
        end

        phibd = factor * philbd + (1.0 - factor) * phirbd
        phibu = factor * philbu + (1.0 - factor) * phirbu

        phifd = factor * philfd + (1.0 - factor) * phirfd
        phifu = factor * philfu + (1.0 - factor) * phirfu

        zbd = factor * zlbd + (1.0 - factor) * zrbd
        zbu = factor * zlbu + (1.0 - factor) * zrbu

        zfd = factor * zlfd + (1.0 - factor) * zrfd
        zfu = factor * zlfu + (1.0 - factor) * zrfu
    end

    # Intepolate in y.
    if sizey == 1
        phid = phibd
        phiu = phibu

        zd = zbd
        zu = zbu
    else
        if yf < yb
            error("Error in interpolate: yf = ", yf, " < yb = ", yb)
        elseif yf == yb
            factor = 0.0
        elseif ylc > yf
            factor = 0.0
        elseif ylc > yb
            factor = (yf - ylc) / (yf - yb)
        else
            factor = 1.0
        end

        phid = factor * phibd + (1.0 - factor) * phifd
        phiu = factor * phibu + (1.0 - factor) * phifu

        zd = factor * zbd + (1.0 - factor) * zfd
        zu = factor * zbu + (1.0 - factor) * zfu
    end

    # Interpolate in z.
    if zu < zd
        error("Error in interpolate: zu = ", zu, " < zd = ", zd)
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    phi = factor * phid + (1.0 - factor) * phiu

    # Return the result.
    return phi
end

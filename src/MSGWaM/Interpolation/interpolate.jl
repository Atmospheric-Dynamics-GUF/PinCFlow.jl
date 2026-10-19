"""
```julia
interpolate(
    input::NamedTuple{
        (:phid, :phiu, :zd, :zu, :zlc),
        <:NTuple{5, <:AbstractFloat},
    },
)::AbstractFloat
```

Perform linear interpolation to `zlc` with the values `(phid, phiu)` at the positions `(zd, zu)`.

The result is given by

```math
\\phi_\\mathrm{C} = f_z \\phi_\\mathrm{D} + (1 - f_z) \\phi_\\mathrm{U}
```

Therein, ``\\mathrm{D}`` and ``\\mathrm{U}`` denote grid points downward and upward of the location of interest (denoted by ``\\mathrm{C}``), respectively, and

```math
f_z = \\begin{cases}
0 & \\mathrm{if} \\quad z_\\mathrm{D} = z_\\mathrm{U},\\\\
1 & \\mathrm{if} \\quad z_\\mathrm{C} < z_\\mathrm{D},\\\\
\\frac{z_\\mathrm{U} - z_\\mathrm{C}}{z_\\mathrm{U} - z_\\mathrm{D}} & \\mathrm{if} \\quad z_\\mathrm{D} \\leq z_\\mathrm{C} \\leq z_\\mathrm{U},\\\\
0 & \\mathrm{if} \\quad z_\\mathrm{U} < z_\\mathrm{C},
\\end{cases}
```

where ``\\left(\\alpha, \\beta, \\gamma\\right) \\in \\left\\{\\left(x, \\mathrm{L}, \\mathrm{R}\\right), \\left(y, \\mathrm{B}, \\mathrm{F}\\right), \\left(z, \\mathrm{D}, \\mathrm{U}\\right)\\right\\}``.

```julia
interpolate(
    state::State,
    input::NamedTuple{
        (
            :philbd,
            :philbu,
            :philfd,
            :philfu,
            :phirbd,
            :phirbu,
            :phirfd,
            :phirfu,
            :zlbd,
            :zlbu,
            :zlfd,
            :zlfu,
            :zrbd,
            :zrbu,
            :zrfd,
            :zrfu,
            :zlc,
            :yb,
            :yf,
            :ylc,
            :xl,
            :xr,
            :xlc,
        ),
        <:NTuple{22, <:AbstractFloat},
    },
)::AbstractFloat
```

Perform trilinear interpolation to `(xlc, ylc, zlc)`, with values from eight surrounding grid points (two zonal positions, two meridional positions and eight vertical positions), and return the result.

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

 1. Interpolation in ``y``:

    ```math
    \\begin{align*}
    \\psi_\\mathrm{D} & = f_y \\psi_\\mathrm{BD} + (1 - f_y) \\psi_\\mathrm{FD},\\\\
    \\psi_\\mathrm{U} & = f_y \\psi_\\mathrm{BU} + (1 - f_y) \\psi_\\mathrm{FU}
    \\end{align*}
    ```

 1. Interpolation in ``z``:

    ```math
    \\phi_\\mathrm{C} = f_z \\phi_\\mathrm{D} + (1 - f_z) \\phi_\\mathrm{U}
    ```

Therein, ``\\mathrm{L}``, ``\\mathrm{R}``, ``\\mathrm{B}``, ``\\mathrm{F}``, ``\\mathrm{D}`` and ``\\mathrm{U}`` denote grid points to the left, to the right, forward, backward, downward and upward of the location of interest (denoted by ``\\mathrm{C}``), respectively, ``\\psi = \\left(\\phi, z\\right)`` and the definitions of ``f_x`` and ``f_y`` are analogous to the above definition of ``f_z``.

# Arguments

  - `state`: Model state.

  - `input`: Values and positions for the interpolation:

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
"""
function interpolate end

function interpolate(
    input::NamedTuple{
        (:phid, :phiu, :zd, :zu, :zlc),
        <:NTuple{5, <:AbstractFloat},
    },
)::AbstractFloat
    (; phid, phiu, zd, zu, zlc) = input

    if zu < zd
        error(
            "Incorrect vertical-interpolation points: zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    return factor * phid + (1.0 - factor) * phiu
end

function interpolate(
    state::State,
    input::NamedTuple{
        (
            :philbd,
            :philbu,
            :philfd,
            :philfu,
            :phirbd,
            :phirbu,
            :phirfd,
            :phirfu,
            :zlbd,
            :zlbu,
            :zlfd,
            :zlfu,
            :zrbd,
            :zrbu,
            :zrfd,
            :zrfu,
            :zlc,
            :yb,
            :yf,
            :ylc,
            :xl,
            :xr,
            :xlc,
        ),
        <:NTuple{22, <:AbstractFloat},
    },
)::AbstractFloat
    (; x_size, y_size) = state.namelists.domain
    (;
        philbd,
        philbu,
        philfd,
        philfu,
        phirbd,
        phirbu,
        phirfd,
        phirfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlc,
        yb,
        yf,
        ylc,
        xl,
        xr,
        xlc,
    ) = input

    # Interpolate in x.
    if x_size == 1
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
            error(
                "Incorrect zonal-interpolation points: xr = ",
                xr,
                " < xl = ",
                xl,
            )
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

    # Interpolate in y.
    if y_size == 1
        phid = phibd
        phiu = phibu

        zd = zbd
        zu = zbu
    else
        if yf < yb
            error(
                "Incorrect meridional-interpolation points: yf = ",
                yf,
                " < yb = ",
                yb,
            )
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

    return interpolate((; phid, phiu, zd, zu, zlc))
end

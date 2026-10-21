"""
```julia
compute_new_grid!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    t::AbstractFloat,
)::Integer
```

Compute the resolved and unresolved topography, Jacobian, metric-tensor elements, physical layers, and minimum layer depth in `grid` and return the number of changes in the resolved topography.

The Jacobian reads

```math
J = \\frac{L_z - h}{L_z} \\frac{\\tilde{z}_{k + 1 / 2} - \\tilde{z}_{k - 1 / 2}}{\\Delta \\hat{z}}
```

and the (non-Cartesian) elements of the metric tensor read

```math
\\begin{align*}
    G^{1 3} & = \\frac{h_{i + 1} - h_{i - 1}}{2 \\Delta \\hat{x}} \\frac{\\tilde{z} - L_z}{L_z - h} \\frac{\\Delta \\hat{z}}{\\tilde{z}_{k + 1 / 2} - \\tilde{z}_{k - 1 / 2}},\\\\
    G^{2 3} & = \\frac{h_{j + 1} - h_{j - 1}}{2 \\Delta \\hat{y}} \\frac{\\tilde{z} - L_z}{L_z - h} \\frac{\\Delta \\hat{z}}{\\tilde{z}_{k + 1 / 2} - \\tilde{z}_{k - 1 / 2}},\\\\
    G^{3 3} & = \\left\\{\\left(\\frac{L_z}{L_z - h}\\right)^2 + \\left(\\frac{\\tilde{z} - L_z}{L_z - h}\\right)^2 \\left[\\left(\\frac{h_{i + 1} - h_{i - 1}}{2 \\Delta \\hat{x}}\\right)^2 + \\left(\\frac{h_{j + 1} - h_{j - 1}}{2 \\Delta \\hat{y}}\\right)^2\\right]\\right\\}\\\\
    & \\quad \\times \\left(\\frac{\\Delta \\hat{z}}{\\tilde{z}_{k + 1 / 2} - \\tilde{z}_{k - 1 / 2}}\\right)^2.
\\end{align*}
```

The vertical layer centers and edges of physical grid are given by

```math
\\begin{align*}
    z_{k + 1 / 2} & = \\frac{L_z - h}{L_z} \\tilde{z}_{k + 1 / 2} + h, & z & = \\frac{L_z - h}{L_z} \\tilde{z} + h.
\\end{align*}
```

This function changes the fields `zcnew` and `zctildenew` of `grid` accordingly but does not change `zc` and `zctilde`.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `t`: Time.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
"""
function compute_new_grid! end

@ivy function compute_new_grid!(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    t::AbstractFloat,
)::Integer
    (; nbz) = namelists.domain
    (; wave_modes, wkb_mode) = namelists.wkb
    (; resolved_topography, unresolved_topography) = namelists.grid
    (; nxx, nyy, nzz, i0, i1, j0, j1, ko, comm) = domain
    (; lref, tref) = constants
    (;
        lz,
        dx,
        dy,
        dz,
        x,
        y,
        zs,
        ztildes,
        hb,
        kh,
        lh,
        hw,
        jac,
        met13,
        met23,
        met33,
        zcnew,
        zctildenew,
        dzcmin,
    ) = grid

    # Compute the resolved topography.
    changes = 0
    for j in j0:j1, i in i0:i1
        hbnew = resolved_topography(x[i] * lref, y[j] * lref, t * tref) / lref

        if hbnew != hb[i, j]
            hb[i, j] = hbnew
            changes += 1
        end
    end
    changes = MPI.Allreduce(changes, +, comm)

    # Compute the unresolved topography.
    if wkb_mode !== :NoWKB
        for j in j0:j1, i in i0:i1, alpha in 1:wave_modes
            (khdim, lhdim, hwdim) =
                unresolved_topography(alpha, x[i] * lref, y[j] * lref, t * tref)
            kh[alpha, i, j] = khdim * lref
            lh[alpha, i, j] = lhdim * lref
            hw[alpha, i, j] = hwdim / lref
        end
    end

    if changes == 0
        return changes
    end

    set_zonal_boundaries_of_field!(hb, namelists, domain)
    set_meridional_boundaries_of_field!(hb, namelists, domain)

    # Set the start index for the computation of the Jacobian and metric tensor.
    kmin = ko == 0 ? 2 : 1
    kmax = nzz

    # Compute the Jacobian.
    for k in kmin:kmax
        jac[:, :, k] .=
            (lz .- hb) ./ lz .* (ztildes[ko + k] .- ztildes[ko + k - 1]) ./ dz
    end
    ko == 0 && (jac[:, :, 1] .= jac[:, :, 2 * nbz])

    # Compute the metric-tensor element 1-3.
    for k in kmin:kmax, j in 1:nyy, i in i0:i1
        met13[i, j, k] =
            (hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx) * (zs[ko + k] - lz) /
            (lz - hb[i, j]) * dz / (ztildes[ko + k] - ztildes[ko + k - 1])
    end
    set_zonal_boundaries_of_field!(met13, namelists, domain)
    ko == 0 && (
        met13[:, :, 1] .=
            met13[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    )

    # Compute the metric-tensor element 2-3.
    for k in 2:nzz, j in j0:j1, i in 1:nxx
        met23[i, j, k] =
            (hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy) * (zs[ko + k] - lz) /
            (lz - hb[i, j]) * dz / (ztildes[ko + k] - ztildes[ko + k - 1])
    end
    set_meridional_boundaries_of_field!(met23, namelists, domain)
    ko == 0 && (
        met23[:, :, 1] .=
            met23[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    )

    # Compute the metric-tensor element 3-3.
    for k in kmin:kmax, j in j0:j1, i in i0:i1
        met33[i, j, k] =
            (
                (lz / (lz - hb[i, j]))^2.0 +
                ((zs[ko + k] - lz) / (lz - hb[i, j]))^2.0 * (
                    ((hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx))^2.0 +
                    ((hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy))^2.0
                )
            ) * (dz / (ztildes[ko + k] - ztildes[ko + k - 1]))^2.0
    end
    ko == 0 && for j in j0:j1, i in i0:i1
        met33[i, j, 1] =
            (
                (lz / (lz - hb[i, j]))^2.0 +
                ((zs[1] - lz) / (lz - hb[i, j]))^2.0 * (
                    ((hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx))^2.0 +
                    ((hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy))^2.0
                )
            ) * (dz / (ztildes[2 * nbz] - ztildes[2 * nbz - 1]))^2.0
    end
    set_zonal_boundaries_of_field!(met33, namelists, domain)
    set_meridional_boundaries_of_field!(met33, namelists, domain)

    # Compute the physical layers.
    for k in 1:nzz
        zctildenew[:, :, k] .= (lz .- hb) ./ lz .* ztildes[ko + k] .+ hb
        zcnew[:, :, k] .= (lz .- hb) ./ lz .* zs[ko + k] .+ hb
    end

    # Compute the minimum layer depth.
    dzcmin[] = mapreduce(
        -,
        min,
        zctildenew[:, :, 2:nzz],
        zctildenew[:, :, 1:(nzz - 1)],
    )
    dzcmin[] = MPI.Allreduce(dzcmin[], min, comm)

    return changes
end

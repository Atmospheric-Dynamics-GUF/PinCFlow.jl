"""
```julia
update_grid!(state::State, time::AbstractFloat)::Nothing
```

Recompute the grid and the atmospheric background fields, interpolate all predictands to the new grid, and shift MS-GWaM's ray volumes such that they are assigned to the correct grid cells.

The transformed vertical wind (``\\hat{w}``) is converted into the Cartesian vertical wind (``w``) before the interpolation, and the result is converted back before the boundary conditions are enforced. In this way, the basis-vector changes are accounted for.

# Arguments

  - `state`: Model state.

  - `time`: Simulation time.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_grid!`](@ref)
"""
function update_grid! end

@ivy function update_grid!(state::State, time::AbstractFloat)::Nothing
    (; namelists, constants, domain, grid, atmosphere) = state
    (; i0, i1, j0, j1, k0, k1, nzz) = domain
    (; zc, zctilde, zcold, zctildeold, dzcmin, dzcminold) = grid
    (; wkb_mode) = namelists.wkb
    (; backups, predictands) = state.variables
    (; tracerbackups, tracerpredictands) = state.tracer
    (; turbulencebackups, turbulencepredictands) = state.turbulence
    (; u, v, w) = predictands
    (; wold) = backups

    # Calculate the Cartesian vertical wind on the old grid.
    dku = 1
    @fuse for k in 1:(nzz - dku), j in j0:j1, i in i0:i1
        wold[i, j, k] = compute_vertical_wind(i, j, k, state)
    end

    @fuse zcold = zc
    @fuse zctildeold = zctilde
    dzcminold[] = dzcmin[]

    changes = compute_grid!(namelists, constants, domain, grid, time)

    changes == 0 && return nothing

    if mapreduce(
        (a, b, c, d) -> max(abs(a - c), abs(b - d)),
        max,
        zc,
        zctilde,
        zcold,
        zctildeold,
    ) > min(dzcmin[], dzcminold[])
        error("Orography is changing too rapidly!")
    end

    compute_background_fields!(namelists, constants, domain, grid, atmosphere)

    # Convert the transformed vertical wind into the Cartesian vertical wind
    # for the interpolation.
    @fuse w[i0:i1, j0:j1, 1:(nzz - dku)] = wold[i0:i1, j0:j1, 1:(nzz - dku)]

    save_backups!(state)

    di = dj = old_grid = true
    for (b, p) in zip(
        (backups, tracerbackups, turbulencebackups),
        (predictands, tracerpredictands, turbulencepredictands),
    )
        for field in fieldnames(typeof(p))
            length(getfield(p, field)) == 0 && continue

            imin = field === :u ? i0 - 1 : i0
            jmin = field === :v ? j0 - 1 : j0
            kmax = (field === :u || field === :v) ? k1 + 1 : k1

            for k in k0:kmax, j in jmin:j1, i in imin:i1
                if field === :u
                    zlc = (zc[i, j, k] + zc[i + 1, j, k]) / 2

                    (kd, ku) = get_levels(i, j, zlc, state; di, old_grid)

                    zd = (zcold[i, j, kd] + zcold[i + 1, j, kd]) / 2
                    zu = (zcold[i, j, ku] + zcold[i + 1, j, ku]) / 2
                elseif field === :v
                    zlc = (zc[i, j, k] + zc[i, j + 1, k]) / 2

                    (kd, ku) = get_levels(i, j, zlc, state; dj, old_grid)

                    zd = (zcold[i, j, kd] + zcold[i, j + 1, kd]) / 2
                    zu = (zcold[i, j, ku] + zcold[i, j + 1, ku]) / 2
                elseif field === :w
                    zlc = zctilde[i, j, k]

                    (kd, ku) = get_half_levels(i, j, zlc, state; dku, old_grid)

                    zd = zctildeold[i, j, kd]
                    zu = zctildeold[i, j, ku]
                else
                    zlc = zc[i, j, k]

                    (kd, ku) = get_levels(i, j, zlc, state; old_grid)

                    zd = zcold[i, j, kd]
                    zu = zcold[i, j, ku]
                end

                phid = getfield(b, Symbol(field, :old))[i, j, kd]
                phiu = getfield(b, Symbol(field, :old))[i, j, ku]

                getfield(p, field)[i, j, k] =
                    interpolate((; phid, phiu, zd, zu, zlc))
            end
        end
    end

    # Convert the Cartesian vertical wind back into the transformed vertical
    # wind.
    @fuse for k in k0:k1, j in j0:j1, i in i0:i1
        uedger = u[i, j, k]
        uuedger = u[i, j, k + 1]
        uedgel = u[i - 1, j, k]
        uuedgel = u[i - 1, j, k + 1]
        vedgef = v[i, j, k]
        vuedgef = v[i, j, k + 1]
        vedgeb = v[i, j - 1, k]
        vuedgeb = v[i, j - 1, k + 1]
        wedgeu = w[i, j, k]

        w[i, j, k] = transform(
            i,
            j,
            k,
            uedger,
            uuedger,
            uedgel,
            uuedgel,
            vedgef,
            vuedgef,
            vedgeb,
            vuedgeb,
            wedgeu,
            Transformed(),
            state,
        )
    end

    set_boundaries!(state, BoundaryPredictands())

    synchronize_compressible_atmosphere!(state, state.variables.predictands)

    if wkb_mode !== :NoWKB
        shift_rays!(state, Val(:SingleColumn))
    end

    nothing
end

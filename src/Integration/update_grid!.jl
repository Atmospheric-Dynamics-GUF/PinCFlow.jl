"""
```julia
update_grid!(state::State, time::AbstractFloat)
```

Recompute the grid and the atmospheric background fields, interpolate all predictands to the new grid, and shift MS-GWaM's ray volumes such that they are assigned to the correct grid cells.

# Arguments

  - `state`: Model state.

  - `time`: Simulation time.
"""
function update_grid! end

function update_grid!(state::State, time::AbstractFloat)
    (; namelists, constants, domain, grid, atmosphere) = state
    (; i0, i1, j0, j1, k0, k1) = domain
    (; zc, zctilde, zcold, zctildeold, dzcmin, dzcminold) = grid
    (; wkb_mode) = namelists.wkb
    (; backups, predictands) = state.variables
    (; tracerbackups, tracerpredictands) = state.tracer
    (; turbulencebackups, turbulencepredictands) = state.turbulence

    zcold .= zc
    zctildeold .= zctilde
    dzcminold[] = dzcmin[]

    changes = compute_new_grid!(namelists, constants, domain, grid, time)

    changes == 0 && return

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

    save_backups!(state)

    for (b, p) in zip(
        (backups, tracerbackups, turbulencebackups),
        (predictands, tracerpredictands, turbulencepredictands),
    )
        for field in fieldnames(typeof(p))
            length(getfield(p, field)) == 0 && continue

            for k in k0:k1, j in j0:j1, i in i0:i1
                if field === :u
                    zlc = (zc[i, j, k] + zc[i + 1, j, k]) / 2

                    (kd, ku) =
                        get_levels(i, j, zlc, state; di = true, old_grid = true)

                    zd = (zcold[i, j, kd] + zcold[i + 1, j, kd]) / 2
                    zu = (zcold[i, j, ku] + zcold[i + 1, j, ku]) / 2
                elseif field === :v
                    zlc = (zc[i, j, k] + zc[i, j + 1, k]) / 2

                    (kd, ku) =
                        get_levels(i, j, zlc, state; dj = true, old_grid = true)

                    zd = (zcold[i, j, kd] + zcold[i, j + 1, kd]) / 2
                    zu = (zcold[i, j, ku] + zcold[i, j + 1, ku]) / 2
                elseif field === :w
                    zlc = zctilde[i, j, k]

                    (kd, ku) =
                        get_half_levels(i, j, zlc, state; old_grid = true)

                    zd = zctildeold[i, j, kd]
                    zu = zctildeold[i, j, ku]
                else
                    zlc = zc[i, j, k]

                    (kd, ku) = get_levels(i, j, zlc, state; old_grid = true)

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

    set_boundaries!(state, BoundaryPredictands())

    synchronize_compressible_atmosphere!(state, state.variables.predictands)

    if wkb_mode !== :NoWKB
        shift_rays!(state, Val(:SingleColumn))
    end

    return
end

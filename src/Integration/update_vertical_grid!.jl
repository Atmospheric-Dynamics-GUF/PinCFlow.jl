"""
```julia
update_vertical_grid!(state::State, time::AbstractFloat)
```

Update the vertical grid by recomputing it, interpolating all predictands to the new grid, recomputing the atmospheric background fields, and shifting MS-GWaM's ray volumes such that they are assigned to the correct grid cells.

# Arguments

  - `state`: Model state.

  - `time`: Simulation time.
"""
function update_vertical_grid! end

function update_vertical_grid!(state::State, time::AbstractFloat)
    (; namelists, constants, domain, grid, atmosphere) = state
    (; i0, i1, j0, j1, k0, k1) = domain
    (; zc, zctilde, zcnew, zctildenew) = grid
    (; wkb_mode) = namelists.wkb
    (; backups, predictands) = state.variables
    (; tracerbackups, tracerpredictands) = state.tracer
    (; turbulencebackups, turbulencepredictands) = state.turbulence

    changes = compute_vertical_grid!(namelists, constants, domain, grid, time)

    if changes == 0
        return
    end

    save_backups!(state)

    for (b, p) in zip(
        (backups, tracerbackups, turbulencebackups),
        (predictands, tracerpredictands, turbulencepredictands),
    )
        for field in fieldnames(typeof(p))
            length(getfield(p, field)) == 0 && continue

            for k in k0:k1, j in j0:j1, i in i0:i1
                if field === :u
                    zlc = (zcnew[i, j, k] + zcnew[i + 1, j, k]) / 2

                    (kd, ku) = get_levels(i, j, zlc, state; di = true)

                    zd = (zc[inew, j, kd] + zc[i + 1, j, kd]) / 2
                    zu = (zc[i, j, ku] + zc[i + 1, j, ku]) / 2
                elseif field === :v
                    zlc = (zcnew[i, j, k] + zcnew[i, j + 1, k]) / 2

                    (kd, ku) = get_levels(i, j, zlc, state; dj = true)

                    zd = (zc[i, j, kd] + zc[i, j + 1, kd]) / 2
                    zu = (zc[i, j, ku] + zc[i, j + 1, ku]) / 2
                elseif field === :w
                    zlc = zctildenew[i, j, k]

                    (kd, ku) = get_half_levels(i, j, zlc, state)

                    zd = zctilde[i, j, kd]
                    zu = zctilde[i, j, ku]
                else
                    zlc = zcnew[i, j, k]

                    (kd, ku) = get_levels(i, j, zlc, state)

                    zd = zc[i, j, kd]
                    zu = zc[i, j, ku]
                end

                phid = getfield(b, Symbol(field, :old))[i, j, kd]
                phiu = getfield(b, Symbol(field, :old))[i, j, ku]

                getfield(p, field)[i, j, k] =
                    interpolate((; phid, phiu, zd, zu, zlc))
            end
        end
    end

    set_boundaries!(state, BoundaryPredictands())

    zc .= zcnew
    zctilde .= zctildenew

    compute_background_fields!(namelists, constants, domain, grid, atmosphere)

    synchronize_compressible_atmosphere!(state, state.variables.predictands)

    if wkb_mode !== :NoWKB
        shift_rays!(state, Val(:SingleColumn))
    end

    return
end

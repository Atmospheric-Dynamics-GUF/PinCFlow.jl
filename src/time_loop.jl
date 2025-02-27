function time_loop!(semi)

    timestep!(semi)

    reconstruction!(semi)

    compute_fluxes!(semi)

end

function timestep!(semi)

    (; equations) = semi
    (; tStepChoice) = equations

    if (tStepChoice == "fix")
        calculate_timestep_fix!(semi)
    elseif (tStepChoice == "cfl")
        calculate_timestep_cfl!(semi)
    end


end

function calculate_timestep_fix!(semi)

    (; dt, tRef, dtMax_dim) = semi.equations

    dt = dtMax_dim / tRef

end

function calculate_timestep_cfl!(semi)

    (; cache, equations, grid) = semi
    (; dt, tRef, dtMax_dim, Re, cfl, lRef, uRef) = equations
    (; var, jac) = cache
    (; nx, ny, nz, dx, dy, dz) = grid

    small = 1.e-20

    u = var.u
    v = var.v
    w = var.w

    uMax = maximum(abs.(u)) + small
    vMax = maximum(abs.(v)) + small
    wMax = maximum(abs.(w)) + small

    dtConv = cfl * min(dx / uMax, dy / vMax, dz / wMax)

    for kz = 1:nz
        for jy = 1:ny
            for ix = 1:nx
                dtConv = min(
                    dtConv,
                    cfl * jac(ix, jy, kz) * dz / (
                        abs(
                            0.5 * (
                                vertWind(ix, jy, kz, semi) + vertWind(ix, jy, kz - 1, semi)
                            ),
                        ) + small
                    ),
                )
            end
        end
    end

    dtVisc = 0.5 * min(dx^2, dy^2, dz^2) * Re

    for kz = 1:nz
        for jy = 1:ny
            for ix = 1:nx
                dtVisc = min(dtVisc, 0.5 * jac(ix, jy, kz)^2 * Re)
            end
        end
    end

    dtMax = dtMax_dim / tRef

    dt = min(dtVisc, dtConv, dtMax)

end




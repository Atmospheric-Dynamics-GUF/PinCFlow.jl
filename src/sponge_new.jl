
function applyUnifiedSponge!(semi, dt)

    applyUnifiedSponge_rho!(semi, dt)
    applyUnifiedSponge_rhop!(semi, dt)
    applyUnifiedSponge_uvw!(semi, dt)

end

function applyUnifiedSponge_rho!(semi, dt)

    (; cache, grid, sponge) = semi
    (; var) = cache
    (; alpha) = sponge
    (; nx, ny, nz) = grid

    rho = var.rho

    for kz = 1:nz
        for jy = 1:ny
            for ix = 1:nx
                alphaSponge = alpha[ix, jy, kz]
                beta = betaSponge(alphaSponge, dt)
                rho[ix, jy, kz] = beta * rho[ix, jy, kz]
            end
        end
    end
end

function applyUnifiedSponge_rhop!(semi, dt)
    # TODO this is exactly the same as for rho, just with rhop
    (; cache, grid, sponge) = semi
    (; var) = cache
    (; alpha) = sponge
    (; nx, ny, nz) = grid

    rhop = var.rhop

    for kz = 1:nz
        for jy = 1:ny
            for ix = 1:nx
                alphaSponge = alpha[ix, jy, kz]
                beta = betaSponge(alphaSponge, dt)
                rhop[ix, jy, kz] = beta * rhop[ix, jy, kz]
            end
        end
    end
end

function applyUnifiedSponge_uvw!(semi, dt)
    # TODO: these are all very similar, only difference 
    # is in alpha!
    applyUnifiedSponge_u!(semi, dt)
    applyUnifiedSponge_v!(semi, dt)
    applyUnifiedSponge_w!(semi, dt)

end

function applyUnifiedSponge_u!(semi, dt)

    (; cache, equations, grid, sponge) = semi
    (; var) = cache
    (; alpha, relax_to_mean) = sponge
    (; nx, ny, nz) = grid
    (; uRef) = equations

    u = var.u
    backgroundFlow_dim = 10
    for kz = 1:nz
        uBG =
            relax_to_mean * sum(u[1:nx, 1:ny, kz]) / (nx * ny) +
            (1.0 - relax_to_mean) * backgroundFlow_dim / uRef * relaxationSponge(semi)

        for jy = 1:ny
            for ix = 1:nx
                alphaSponge = alphaSponge_x(ix, jy, kz, alpha)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, u[ix, jy, kz], uBG)
            end
        end
    end
end

function applyUnifiedSponge_v!(semi, dt)

    (; cache, equations, grid, sponge) = semi
    (; var) = cache
    (; alpha, relax_to_mean) = sponge
    (; nx, ny, nz) = grid
    (; uRef) = equations

    v = var.v
    backgroundFlow_dim = 0.0
    for kz = 1:nz
        vBG =
            relax_to_mean * sum(v[1:nx, 1:ny, kz]) / (nx * ny) +
            (1.0 - relax_to_mean) * backgroundFlow_dim / uRef * relaxationSponge(semi)

        for jy = 1:ny
            for ix = 1:nx
                alphaSponge = alphaSponge_y(ix, jy, kz, alpha)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, v[ix, jy, kz], vBG)
            end
        end
    end
end

function applyUnifiedSponge_w!(semi, dt)

    (; cache, equations, grid, sponge) = semi
    (; var, jac) = cache
    (; alpha, relax_to_mean) = sponge
    (; nx, ny, nz) = grid
    (; uRef) = equations
    backgroundFlow_dim = 0.0

    w = var.w

    for kz = 1:nz
        wBG =
            relax_to_mean * sum(w[1:nx, 1:ny, kz]) / (nx * ny) +
            (1.0 - relax_to_mean) * backgroundFlow_dim / uRef * relaxationSponge(semi)

        for jy = 1:ny
            for ix = 1:nx
                alphaSponge = alphaSponge_z(ix, jy, kz, alpha, jac)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, w[ix, jy, kz], wBG)
            end
        end
    end
end

function alphaSponge_x(ix, jy, kz, alpha)

    alphaSponge = 0.5 * (alpha[ix, jy, kz] + alpha[ix+1, jy, kz])

    return alphaSponge
end

function alphaSponge_y(ix, jy, kz, alpha)

    alphaSponge = 0.5 * (alpha[ix, jy, kz] + alpha[ix, jy+1, kz])

    return alphaSponge
end

function alphaSponge_z(ix, jy, kz, alpha, jac)

    alphaSponge =
        (jac(ix, jy, kz + 1) * alpha[ix, jy, kz] + jac(ix, jy, kz) * alpha[ix, jy, kz+1]) /
        (jac(ix, jy, kz) + jac(ix, jy, kz + 1))

    return alphaSponge
end

function betaSponge(alphaSponge, dt)

    return 1.0 / (1.0 + alphaSponge * dt)

end

function unifiedSponge_new!(beta, uvw, uvwBG)

    uvw = (1.0 - beta) * uvwBG + beta * uvw

end

function relaxationSponge(semi)

    (; sponge, equations) = semi
    (; relaxation_amplitude, relaxation_period) = sponge
    (; tRef) = equations # TODO needs time

    time = 0.0

    if (relaxation_period > 0.0)
        return (
            1.0 + relaxation_amplitude * sin(2.0 * pi * time / relaxation_period * tRef)
        )
    else
        return 1.0
    end

end

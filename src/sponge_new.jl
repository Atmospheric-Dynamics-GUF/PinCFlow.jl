
function applyUnifiedSponge!(semi, dt)
    applyUnifiedSponge_rho!(semi, dt)
    applyUnifiedSponge_rhop!(semi, dt)
    applyUnifiedSponge_uvw!(semi, dt)
end

function applyUnifiedSponge_rho!(model, dt)
    @trixi_timeit timer() "Sponge layer rho" begin
    #! format: noindent
    (; nx, ny, nz) = model.domain
    alpha = model.atmosphere.sponge.alpha
    rho = model.variables.prognostic_fields.rho

    for kz in 1:nz
        for jy in 1:ny
            for ix in 1:nx
                alphaSponge = alpha[ix, jy, kz]
                beta = betaSponge(alphaSponge, dt)
                rho[ix, jy, kz] = beta * rho[ix, jy, kz]
            end
        end
    end
    end # timer
end

function applyUnifiedSponge_rhop!(model, dt)
    @trixi_timeit timer() "Sponge layer rhop" begin
    #! format: noindent
    # TODO this is exactly the same as for rho, just with rhop
    (; nx, ny, nz) = model.domain
    alpha = model.atmosphere.sponge.alpha
    rhop = model.variables.prognostic_fields.rho

    for kz in 1:nz
        for jy in 1:ny
            for ix in 1:nx
                alphaSponge = alpha[ix, jy, kz]
                beta = betaSponge(alphaSponge, dt)
                rhop[ix, jy, kz] = beta * rhop[ix, jy, kz]
            end
        end
    end
    end # timer
end

function applyUnifiedSponge_uvw!(semi, dt)
    @trixi_timeit timer() "sponge_uvw" begin
    #! format: noindent
    # TODO: these are all very similar, only difference
    # is in alpha!
    applyUnifiedSponge_u!(semi, dt)
    applyUnifiedSponge_v!(semi, dt)
    applyUnifiedSponge_w!(semi, dt)
    end # timer
end

function applyUnifiedSponge_u!(model, dt)
    (; uref) = model.constants
    (; nx, ny, nz) = model.domain
    (; alpha, relax_to_mean) = model.atmosphere.sponge
    u = model.variables.prognostic_fields.u
    # TODO:
    backgroundFlow_dim = 10
    for kz in 1:nz
        uBG = relax_to_mean * sum(u[1:nx, 1:ny, kz]) / (nx * ny) +
              (1.0 - relax_to_mean) * backgroundFlow_dim / uref * relaxationSponge(model)

        for jy in 1:ny
            for ix in 1:nx
                alphaSponge = alphaSponge_x(ix, jy, kz, alpha)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, u, uBG, ix, jy, kz)
            end
        end
    end
end

function applyUnifiedSponge_v!(model, dt)

    (; uref) = model.constants
    (; nx, ny, nz) = model.domain
    (; alpha, relax_to_mean) = model.atmosphere.sponge
    v = model.variables.prognostic_fields.v
    # TODO:
    backgroundFlow_dim = 0.0
    for kz in 1:nz
        vBG = relax_to_mean * sum(v[1:nx, 1:ny, kz]) / (nx * ny) +
              (1.0 - relax_to_mean) * backgroundFlow_dim / uref * relaxationSponge(model)

        for jy in 1:ny
            for ix in 1:nx
                alphaSponge = alphaSponge_y(ix, jy, kz, alpha)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, v, vBG, ix, jy, kz)
            end
        end
    end
end

function applyUnifiedSponge_w!(model, dt)
    backgroundFlow_dim = 0.0
    (; uref) = model.constants
    (; nx, ny, nz) = model.domain
    (; alpha, relax_to_mean) = model.atmosphere.sponge
    w = model.variables.prognostic_fields.w

    jac = model.grid.jac

    for kz in 1:nz
        wBG = relax_to_mean * sum(w[1:nx, 1:ny, kz]) / (nx * ny) +
              (1.0 - relax_to_mean) * backgroundFlow_dim / uref * relaxationSponge(model)

        for jy in 1:ny
            for ix in 1:nx
                alphaSponge = alphaSponge_z(ix, jy, kz, alpha, jac)
                beta = betaSponge(alphaSponge, dt)
                unifiedSponge_new!(beta, w, wBG, ix, jy, kz)
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
    alphaSponge = (jac(ix, jy, kz + 1) * alpha[ix, jy, kz] +
                   jac(ix, jy, kz) * alpha[ix, jy, kz+1]) /
                  (jac(ix, jy, kz) + jac(ix, jy, kz + 1))

    return alphaSponge
end

function betaSponge(alphaSponge, dt)
    return 1.0 / (1.0 + alphaSponge * dt)
end

function unifiedSponge_new!(beta, uvw, uvwBG, indices...)
    uvw[indices...] = (1.0 - beta) * uvwBG + beta * uvw[indices...]
end

function relaxationSponge(model)

    (; relaxation_amplitude, relaxation_period) = model.atmosphere.sponge
    tref = model.constants.tref

    time = 0.0

    if (relaxation_period > 0.0)
        return (1.0 +
                relaxation_amplitude * sin(2.0 * pi * time / relaxation_period * tref))
    else
        return 1.0
    end
end

function time_loop!(semi)

    timestep!(semi)

    reconstruction!(semi)

    calculate_fluxes!(semi)

    rk3!(semi)

end

function timestep!(semi)

    (; tStepChoice) = semi 

    calculate_timestep!(semi, tStepChoice)

end

function calculate_timestep!(semi)

    (; dt, tRef, dtMax_dim) = semi 

    dt = dtMax_dim / tRef

end

function calculate_timestep!(semi)

    (; dt, tRef, dtMax_dim, cache, Re, grid) = semi 
    (; var) = cache
    (; nx, ny, nz, dx, dy, dz) = grid

    small = 1.e-20

    u =  var.u 
    v =  var.v 
    w =  var.w

    uMax = maximum(abs(u)) + small 
    vMax = maximum(abs(v)) + small
    wMax = maximum(abs(w)) + small

    dtConv = cfl * min(dx / uMax, dy / vMax, dz / wMax)

    for kz in 1 : nz 
        for jy in 1 : ny 
            for ix in 1 : nx 
                dtConv = min(dtConv, cfl * jac(ix, jy, kz) * dz / 
                    (abs(0.5 * (vertWindTFC(ix, jy, kz) + vertWindTFC(ix, jy, kz - 1)))
                    + small))
            end
        end
    end

    dtVisc = 0.5 * min(dx ^ 2, dy ^ 2, dz ^ 2) / Re

    for kz in 1 : nz 
        for jy in 1 : ny 
            for ix in 1 : nx 
                dtVisc = min(dtVisc, 0.5 * jac(ix, jy, kz) ^ 2 * Re)
            end
        end
    end

    dtMax = dtMax_dim / tRef

    dt = min(dtVisc, dtConv, dtMax)
    
end



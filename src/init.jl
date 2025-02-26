function initialize_variables!(semi)

    (; cache, grid, refval) = semi
    (; x, y, z, zTFC, nx, ny, nz, nbx, nby, nbz) = grid 
    (; var, rhoStrat, pStrat, thetaStrat, bvsStrat) = cache 
    (; lRef, uRef) = refval

    u = @view var.u 
    v = @view var.v
    w = @view var.w
    rho = @view var.rho
    exner = @view var.exner
    
    u0 = 10. # this should be a namelist parameter
    v0 = 0.
    w0 = 0.

    u .= u0 / uRef
    v .= v0 / uRef
    w .= w0 / uRef

    rho .= 0.0
    exner .= 0.0

end

function initialize!(semi)

    initialize_atmosphere!(semi)

    initialize_variables!(semi)

end
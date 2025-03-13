function setBoundary!(model::Model)
    @trixi_timeit timer() "Set boundary" begin
    #! format: noindent
    setBoundary_x!(model, model.grid.xboundary)
    setBoundary_y!(model, model.grid.yboundary)
    setBoundary_z!(model, model.grid.zboundary)
    end
end

function setBoundary_flux!(model)
    nz = model.domain.nz
    flux = model.fluxes

    flux.rho[:, :, 0, 3] .= 0.0
    flux.rho[:, :, nz, 3] .= 0.0
    flux.rhop[:, :, 0, 3] .= 0.0
    flux.rhop[:, :, nz, 3] .= 0.0

    flux.u[:, :, 0, 3] .= 0.0
    flux.u[:, :, nz, 3] .= 0.0
    flux.v[:, :, 0, 3] .= 0.0
    flux.v[:, :, nz, 3] .= 0.0
    flux.w[:, :, -1, 3] .= 0.0
    flux.w[:, :, nz, 3] .= 0.0
end

function setBoundary_x!(model, boundary::PeriodicBC)
    (; rho, rhop, u, v, w, pip) = model.variables.prognostic_fields
    (; nx, nbx) = model.domain

    for i in 1:nbx
        set_periodic_value_cell_x!(rho, i, nx)
        set_periodic_value_cell_x!(rhop, i, nx)
    end

    @views u[0, :, :] .= u[nx, :, :]

    for i in 1:nbx
        set_periodic_value_edge_x!(u, i, nx)
        set_periodic_value_cell_x!(v, i, nx)
        set_periodic_value_cell_x!(w, i, nx)
    end

    # TODO: why are we not looping over x indices as for rho and rhop?
    set_periodic_value_cell_x!(pip, 1, nx)
    ## Missing RhoTilde, but it should be removed according to PincFlow.f90
end

function set_periodic_value_cell_x!(var, i, nx)
    @views var[nx+i, :, :] .= var[i, :, :]
    @views var[-i+1, :, :] .= var[nx-i+1, :, :]
end

function set_periodic_value_edge_x!(var, i, nx)
    @views var[nx+i, :, :] .= var[i, :, :]
    @views var[-i, :, :] .= var[nx-i, :, :]
end

function setBoundary_y!(model, boundary::PeriodicBC)
    (; rho, rhop, u, v, w, pip) = model.variables.prognostic_fields
    (; ny, nby) = model.domain

    for j in 1:nby
        set_periodic_value_cell_y!(rho, j, ny)
        set_periodic_value_cell_y!(rhop, j, ny)
    end

    @views v[:, 0, :] .= v[:, ny, :]

    for j in 1:nby
        set_periodic_value_cell_y!(u, j, ny)
        set_periodic_value_edge_y!(v, j, ny)
        set_periodic_value_cell_y!(w, j, ny)
    end

    set_periodic_value_cell_y!(pip, 1, ny)
end

function set_periodic_value_cell_y!(var, j, ny)
    @views var[:, ny+j, :] .= var[:, j, :]
    @views var[:, -j+1, :] .= var[:, ny-j+1, :]
end

function set_periodic_value_edge_y!(var, j, ny)
    @views var[:, ny+j, :] .= var[:, j, :]
    @views var[:, -j, :] .= var[:, ny-j, :]
end

function setBoundary_z!(model, boundary::PeriodicBC)
    println("in periodic")
    (; rho, rhop, u, v, w, pip) = model.variables.prognostic_fields
    (; nz, nbz) = model.domain

    # TODO: This should be improved...

    for k in 1:nbz
        set_periodic_value_cell_z!(rho, k, nz)
        set_periodic_value_cell_z!(rhop, k, nz)
    end

    @views w[:, :, 0] .= w[:, :, nz]

    for k in 1:nbz
        set_periodic_value_cell_z!(u, k, nz)
        set_periodic_value_cell_z!(v, k, nz)
        set_periodic_value_edge_z!(w, k, nz)
    end

    set_periodic_value_cell_z!(pip, 1, nz)
end

function set_periodic_value_cell_z!(var, k, nz)
    @views var[:, :, nz+k] .= var[:, :, k]
    @views var[:, :, -k+1] .= var[:, :, nz-k+1]
end

function set_periodic_value_edge_z!(var, k, nz)
    @views var[:, :, nz+k] .= var[:, :, k]
    @views var[:, :, -k] .= var[:, :, nz-k]
end

function setBoundary_z!(model, boundary::SolidWallBC)
    println("in solid")
    (; rho, rhop, u, v, w, pip) = model.variables.prognostic_fields
    (; nz, nbz) = model.domain

    ## This should be removed!!
    for k in 1:nbz
        set_reflective_value_cell_z!(rho, k, nz)
        set_reflective_value_cell_z!(rhop, k, nz)
    end

    w[:, :, 0] .= 0.0
    w[:, :, nz] .= 0.0

    for k in 1:nbz
        set_reflective_value_edge_z!(w, k, nz)
        set_non_reflective_value_cell_z!(u, k, nz)
        set_non_reflective_value_cell_z!(v, k, nz)
    end

    set_non_reflective_value_cell_z!(pip, 1, nz)
end

function set_reflective_value_cell_z!(var, k, nz)
    for i in axes(var, 1), j in axes(var, 2)
        var[i, j, nz+k] = -var[i, j, nz-k+1]
        var[i, j, -k+1] = -var[i, j, k]
    end
end

function set_reflective_value_edge_z!(var, k, nz)
    for i in axes(var, 1), j in axes(var, 2)
        var[i, j, nz+k] = -var[i, j, nz-k]
        var[i, j, -k] = -var[i, j, k]
    end
end

function set_non_reflective_value_cell_z!(var, k, nz)
    @views var[:, :, nz+k] .= var[:, :, nz-k+1]
    @views var[:, :, -k+1] .= var[:, :, k]
end

function set_non_reflective_value_edge_z!(var, k, nz)
    @views var[:, :, nz+k] .= var[:, :, nz-k]
    @views var[:, :, -k] .= var[:, :, k]
end

function setBoundary!(semi)
    setBoundary_x!(semi, semi.boundary_conditions.boundary_x)
    setBoundary_y!(semi, semi.boundary_conditions.boundary_y)
    setBoundary_z!(semi, semi.boundary_conditions.boundary_z)
end

function setBoundary_flux!(semi)
    (; cache, grid) = semi
    (; nz) = grid
    (; flux) = cache

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

function setBoundary_x!(semi, boundary::PeriodicBC)
    (; cache, grid) = semi
    (; rho, rhop, u, v, w, exner) = cache.var
    (; nx, nbx) = grid

    for i in 1:nbx
        set_periodic_value_cell_x!(rho, i, nx)
        set_periodic_value_cell_x!(rhop, i, nx)
    end

    u[0, :, :] .= u[nx, :, :]

    for i in 1:nbx
        set_periodic_value_edge_x!(u, i, nx)
        set_periodic_value_cell_x!(v, i, nx)
        set_periodic_value_cell_x!(w, i, nx)
    end

    set_periodic_value_cell_x!(exner, 1, nx)
    ## Missing RhoTilde, but it should be removed according to PincFlow.f90
end

function set_periodic_value_cell_x!(var, i, nx)
    var[nx + i, :, :] .= var[i, :, :]
    var[-i + 1, :, :] .= var[nx - i + 1, :, :]
end

function set_periodic_value_edge_x!(var, i, nx)
    var[nx + i, :, :] .= var[i, :, :]
    var[-i, :, :] .= var[nx - i, :, :]
end

function setBoundary_y!(semi, boundary::PeriodicBC)
    (; cache, grid) = semi
    (; rho, rhop, u, v, w, exner) = cache.var
    (; ny, nby) = grid

    for j in 1:nby
        set_periodic_value_cell_y!(rho, j, ny)
        set_periodic_value_cell_y!(rhop, j, ny)
    end

    v[:, 0, :] .= v[:, ny, :]

    for j in 1:nby
        set_periodic_value_cell_y!(u, j, ny)
        set_periodic_value_edge_y!(v, j, ny)
        set_periodic_value_cell_y!(w, j, ny)
    end

    set_periodic_value_cell_y!(exner, 1, ny)
end

function set_periodic_value_cell_y!(var, j, ny)
    var[:, ny + j, :] .= var[:, j, :]
    var[:, -j + 1, :] .= var[:, ny - j + 1, :]
end

function set_periodic_value_edge_y!(var, j, ny)
    var[:, ny + j, :] .= var[:, j, :]
    var[:, -j, :] .= var[:, ny - j, :]
end

function setBoundary_z!(semi, boundary::PeriodicBC)
    (; cache, grid) = semi
    (; rho, rhop, u, v, w, exner) = cache.var
    (; nz, nbz) = grid

    for k in 1:nbz
        set_periodic_value_cell_z!(rho, k, nz)
        set_periodic_value_cell_z!(rhop, k, nz)
    end

    w[:, :, 0] .= v[:, :, nz]

    for k in 1:nbz
        set_periodic_value_cell_z!(u, k, nz)
        set_periodic_value_cell_z!(v, k, nz)
        set_periodic_value_edge_z!(w, k, nz)
    end

    set_periodic_value_cell_z!(exner, 1, nz)
end

function set_periodic_value_cell_z!(var, k, nz)
    var[:, :, nz + k] .= var[:, :, k]
    var[:, :, -k + 1] .= var[:, :, nz - k + 1]
end

function set_periodic_value_edge_z!(var, k, nz)
    var[:, :, nz + k] .= var[:, :, k]
    var[:, -k, :] .= var[:, :, nz - k]
end

function setBoundary_z!(semi, boundary::SolidWallBC)
    (; cache, grid) = semi
    (; rho, rhop, u, v, w, exner) = cache.var
    (; nz, nbz) = grid

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

    set_non_reflective_value_cell_z!(exner, 1, nz)
end

function set_reflective_value_cell_z!(var, k, nz)
    var[:, :, nz + k] .= -var[:, :, nz - k + 1]
    var[:, :, -k + 1] .= -var[:, :, k]
end

function set_reflective_value_edge_z!(var, k, nz)
    var[:, :, nz + k] .= -var[:, :, nz - k]
    var[:, :, -k] .= -var[:, :, k]
end

function set_non_reflective_value_cell_z!(var, k, nz)
    var[:, :, nz + k] .= var[:, :, nz - k + 1]
    var[:, :, -k + 1] .= var[:, :, k]
end

function set_non_reflective_value_edge_z!(var, k, nz)
    var[:, :, nz + k] .= var[:, :, nz - k]
    var[:, :, -k] .= var[:, :, k]
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    boundaries::SolidWallBoundaries,
)

    # Get all necessary fields.
    (; nbz) = state.namelists.domain
    (; k0, k1) = state.domain
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    # Set the density fluctuations at the boundaries to zero.
    for k in 1:nbz
        @views rho[:, :, k0 - k] .= -rho[:, :, k0 + k - 1]
        @views rho[:, :, k1 + k] .= -rho[:, :, k1 - k + 1]
    end

    # Set the density fluctuations at the boundaries to zero.
    for k in 1:nbz
        @views rhop[:, :, k0 - k] .= -rhop[:, :, k0 + k - 1]
        @views rhop[:, :, k1 + k] .= -rhop[:, :, k1 - k + 1]
    end

    # Set the vertical wind at the boundaries to zero.
    w[:, :, k0 - 1] .= 0.0
    w[:, :, k1] .= 0.0
    for k in 1:nbz
        @views w[:, :, k0 - k] .= -w[:, :, k0 + k - 2]
        @views w[:, :, k1 + k] .= -w[:, :, k1 - k]
    end

    # Set the horizontal-wind gradient at the boundaries to zero.
    for k in 1:nbz
        @views u[:, :, k0 - k] .= u[:, :, k0 + k - 1]
        @views u[:, :, k1 + k] .= u[:, :, k1 - k + 1]
        @views v[:, :, k0 - k] .= v[:, :, k0 + k - 1]
        @views v[:, :, k1 + k] .= v[:, :, k1 - k + 1]
    end

    # Set the pressure gradient at the boundaries to zero.
    @views pip[:, :, k0 - 1] .= pip[:, :, k0]
    @views pip[:, :, k1 + 1] .= pip[:, :, k1]

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    boundaries::SolidWallBoundaries,
)

    # Get all necessary fields.
    (; k0, k1) = state.domain
    (; phirho, phirhop, phiu, phiv, phiw) = state.variables.fluxes

    # Set all vertical boundary fluxes to zero.

    # Density
    phirho[:, :, k0 - 1, 3] .= 0.0
    phirho[:, :, k1, 3] .= 0.0

    # Density fluctuations
    phirhop[:, :, k0 - 1, 3] .= 0.0
    phirhop[:, :, k1, 3] .= 0.0

    # Zonal momentum
    phiu[:, :, k0 - 1, 3] .= 0.0
    phiu[:, :, k1, 3] .= 0.0

    # Meridional momentum
    phiv[:, :, k0 - 1, 3] .= 0.0
    phiv[:, :, k1, 3] .= 0.0

    # Vertical momentum
    phiw[:, :, k0 - 2, 3] .= 0.0
    phiw[:, :, k1, 3] .= 0.0

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    boundaries::AbstractBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, boundaries, wkb_mode)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    boundaries::SolidWallBoundaries,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; uw, vw, e) = state.wkb.integrals

    for jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        uw[ix, jy, k0 - 1] = uw[ix, jy, k0]
        uw[ix, jy, k1 + 1] = uw[ix, jy, k1]
        vw[ix, jy, k0 - 1] = vw[ix, jy, k0]
        vw[ix, jy, k1 + 1] = vw[ix, jy, k1]
        e[ix, jy, k0 - 1] = e[ix, jy, k0]
        e[ix, jy, k1 + 1] = e[ix, jy, k1]
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    boundaries::SolidWallBoundaries,
    wkb_mode::MultiColumn,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; uu, uv, uw, vv, vw, etx, ety, utheta, vtheta, e) = state.wkb.integrals

    for jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        uu[ix, jy, k0 - 1] = uu[ix, jy, k0]
        uu[ix, jy, k1 + 1] = uu[ix, jy, k1]
        uv[ix, jy, k0 - 1] = uv[ix, jy, k0]
        uv[ix, jy, k1 + 1] = uv[ix, jy, k1]
        uw[ix, jy, k0 - 1] = uw[ix, jy, k0]
        uw[ix, jy, k1 + 1] = uw[ix, jy, k1]
        vv[ix, jy, k0 - 1] = vv[ix, jy, k0]
        vv[ix, jy, k1 + 1] = vv[ix, jy, k1]
        vw[ix, jy, k0 - 1] = vw[ix, jy, k0]
        vw[ix, jy, k1 + 1] = vw[ix, jy, k1]
        etx[ix, jy, k0 - 1] = etx[ix, jy, k0]
        etx[ix, jy, k1 + 1] = etx[ix, jy, k1]
        ety[ix, jy, k0 - 1] = ety[ix, jy, k0]
        ety[ix, jy, k1 + 1] = ety[ix, jy, k1]
        utheta[ix, jy, k0 - 1] = utheta[ix, jy, k0]
        utheta[ix, jy, k1 + 1] = utheta[ix, jy, k1]
        vtheta[ix, jy, k0 - 1] = vtheta[ix, jy, k0]
        vtheta[ix, jy, k1 + 1] = vtheta[ix, jy, k1]
        e[ix, jy, k0 - 1] = e[ix, jy, k0]
        e[ix, jy, k1 + 1] = e[ix, jy, k1]
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    boundaries::AbstractBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, boundaries, wkb_mode)
    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    boundaries::SolidWallBoundaries,
    wkb_mode::Union{SteadyState, SingleColumn},
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dudt, dvdt) = state.wkb.tendencies

    for jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        dudt[ix, jy, k0 - 1] = dudt[ix, jy, k0]
        dudt[ix, jy, k1 + 1] = dudt[ix, jy, k1]
        dvdt[ix, jy, k0 - 1] = dvdt[ix, jy, k0]
        dvdt[ix, jy, k1 + 1] = dvdt[ix, jy, k1]
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    boundaries::SolidWallBoundaries,
    wkb_mode::MultiColumn,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    for jy in (j0 - 1):(j1 + 1), ix in (i0 - 1):(i1 + 1)
        dudt[ix, jy, k0 - 1] = dudt[ix, jy, k0]
        dudt[ix, jy, k1 + 1] = dudt[ix, jy, k1]
        dvdt[ix, jy, k0 - 1] = dvdt[ix, jy, k0]
        dvdt[ix, jy, k1 + 1] = dvdt[ix, jy, k1]
        dthetadt[ix, jy, k0 - 1] = dthetadt[ix, jy, k0]
        dthetadt[ix, jy, k1 + 1] = dthetadt[ix, jy, k1]
    end

    return
end

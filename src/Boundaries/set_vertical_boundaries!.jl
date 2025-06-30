"""
    set_vertical_boundaries!(state::State, variables::BoundaryPredictands, zboundaries::SolidWallBoundaries)

Enforce vertical boundary conditions for all predictand fields.

The symmetry conditions are as follows:
- Density-fluctuation fields (`rho`, `rhop`): point reflection (-)
- Vertical velocity (`w`): point reflection (-) on the staggered grid
- Horizontal velocities (`u`, `v`): line reflection (+)
- Exner-pressure fluctuations (`pip`): line reflection (+)
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    zboundaries::SolidWallBoundaries,
)
    (; namelists, domain) = state
    (; model, zboundaries) = namelists.setting
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    set_vertical_boundaries_of_field!(rho, namelists, domain, zboundaries, -)
    set_vertical_boundaries_of_field!(rhop, namelists, domain, zboundaries, -)

    set_vertical_boundaries_of_field!(
        w,
        namelists,
        domain,
        zboundaries,
        -;
        staggered = true,
    )

    set_vertical_boundaries_of_field!(u, namelists, domain, zboundaries, +)
    set_vertical_boundaries_of_field!(v, namelists, domain, zboundaries, +)

    set_vertical_boundaries_of_field!(pip, namelists, domain, zboundaries, +)

    set_compressible_vertical_boundaries!(state, variables, model)

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryReconstructions, zboundaries::SolidWallBoundaries)

Enforce vertical boundary conditions for all reconstruction fields.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    zboundaries::SolidWallBoundaries,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; reconstructions) = state.variables

    for field in fieldnames(Reconstructions)
        set_vertical_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
            zboundaries,
        )
    end

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryFluxes, zboundaries::SolidWallBoundaries)

Set the vertical fluxes at the vertical boundaries to zero (in `SolidWallBoundaries` configurations).
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    zboundaries::SolidWallBoundaries,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables
    (; model) = state.namelists.setting

    # Set all vertical boundary fluxes to zero.

    if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    if ko + nzz == sizezz
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    set_compressible_vertical_boundaries!(state, variables, model)

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWIntegrals, zboundaries::SolidWallBoundaries)

Enforce vertical boundary conditions for gravity-wave-integral fields, dispatching based on WKB mode.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    zboundaries::SolidWallBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, zboundaries, wkb_mode)
    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWIntegrals, zboundaries::SolidWallBoundaries, wkb_mode::AbstractWKBMode)

Enforce vertical boundary conditions for gravity-wave-integral fields needed in `SingleColumn` and `SteadyState` configurations, using line reflection.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    zboundaries::SolidWallBoundaries,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_vertical_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain,
            zboundaries,
            +;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWIntegrals, zboundaries::SolidWallBoundaries, wkb_mode::MultiColumn)

Enforce vertical boundary conditions for gravity-wave-integral fields needed in `MultiColumn` configurations, using line reflection.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    zboundaries::SolidWallBoundaries,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_vertical_boundaries_of_field!(
            getfield(integrals, field),
            namelists,
            domain,
            zboundaries,
            +;
            layers = (1, 1, 1),
        )
    end

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWTendencies, zboundaries::SolidWallBoundaries)

Enforce vertical boundary conditions for gravity-wave-tendency fields, dispatching based on WKB mode.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    zboundaries::SolidWallBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    set_vertical_boundaries!(state, variables, zboundaries, wkb_mode)
    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWTendencies, zboundaries::SolidWallBoundaries, wkb_mode::AbstractWKBMode)

Enforce vertical boundary conditions for gravity-wave-tendency fields needed in `SingleColumn` and `SteadyState` configurations, using line reflection.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    zboundaries::SolidWallBoundaries,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt)
        set_vertical_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
            zboundaries,
            +,
        )
    end

    return
end

"""
    set_vertical_boundaries!(state::State, variables::BoundaryGWTendencies, zboundaries::SolidWallBoundaries, wkb_mode::MultiColumn)

Enforce vertical boundary conditions for gravity-wave-tendency fields needed in `MultiColumn` configurations, using line reflection.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    zboundaries::SolidWallBoundaries,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; tendencies) = state.wkb

    for field in (:dudt, :dvdt, :dthetadt)
        set_vertical_boundaries_of_field!(
            getfield(tendencies, field),
            namelists,
            domain,
            zboundaries,
            +,
        )
    end

    return
end

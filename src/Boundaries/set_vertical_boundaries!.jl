function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    boundaries::AbstractBoundaries,
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

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    boundaries::AbstractBoundaries,
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

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    boundaries::SolidWallBoundaries,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables
    (; model) = state.namelists.setting

    # Set all vertical boundary fluxes to zero.

    if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    set_compressible_vertical_boundaries!(state, variables, model)

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
    boundaries::AbstractBoundaries,
    wkb_mode::AbstractWKBMode,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; integrals) = state.wkb

    for field in (:uw, :vw, :e)
        set_vertical_boundaries_of_reduced_field!(
            getfield(integrals, field),
            namelists,
            domain,
            zboundaries,
            +,
        )
    end

    return
end

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    boundaries::AbstractBoundaries,
    wkb_mode::MultiColumn,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; integrals) = state.wkb

    for field in (:uu, :uv, :uw, :vv, :vw, :etx, :ety, :utheta, :vtheta, :e)
        set_vertical_boundaries_of_reduced_field!(
            getfield(integrals, field),
            namelists,
            domain,
            zboundaries,
            +,
        )
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
    boundaries::AbstractBoundaries,
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

function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    boundaries::AbstractBoundaries,
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

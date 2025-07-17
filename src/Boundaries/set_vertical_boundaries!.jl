"""
    set_vertical_boundaries!(state, variables::BoundaryPredictands, boundaries)

Set vertical boundaries for predictand fields with appropriate symmetry conditions:

  - Density fields (rho, rhop): antisymmetric (-)
  - Vertical velocity (w): antisymmetric (-), staggered grid handling
  - Horizontal velocities (u, v): symmetric (+)
  - Pressure perturbation (pip): symmetric (+)
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    boundaries::AbstractBoundaries,
)
    (; namelists, domain) = state
    (; model, zboundaries) = namelists.setting
    (; rho, rhop, u, v, w, pip) = state.variables.predictands
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

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

    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

"""
    set_vertical_boundaries!(state, variables::BoundaryReconstructions, boundaries)

Set vertical boundaries for all reconstruction fields.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    boundaries::AbstractBoundaries,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; reconstructions) = state.variables
    (; tracersetup) = namelists.tracer
    (; icesetup) = namelists.ice
    (; turbulencesetup) = namelists.turbulence

    for field in fieldnames(Reconstructions)
        set_vertical_boundaries_of_field!(
            getfield(reconstructions, field),
            namelists,
            domain,
            zboundaries,
        )
    end

    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

"""
    set_vertical_boundaries!(state, variables::BoundaryFluxes, boundaries::SolidWallBoundaries)

Set vertical flux boundaries to zero at solid walls (top and bottom domain boundaries).
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    boundaries::SolidWallBoundaries,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; fluxes) = state.variables
    (; model) = state.namelists.setting
    (; tracersetup) = state.namelists.tracer
    (; icesetup) = state.namelists.ice
    (; turbulencesetup) = state.namelists.turbulence

    # Set all vertical boundary fluxes to zero.

    if ko == 0
        for field in (:phirho, :phirhop, :phiu, :phiv, :phitheta)
            getfield(fluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
        fluxes.phiw[:, :, k0 - 2, 3] .= 0.0
    end

    if ko + nzz == sizezz
        for field in (:phirho, :phirhop, :phiu, :phiv, :phiw, :phitheta)
            getfield(fluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    set_compressible_vertical_boundaries!(state, variables, model)
    set_tracer_vertical_boundaries!(state, variables, tracersetup)
    set_ice_vertical_boundaries!(state, variables, icesetup)
    set_turbulence_vertical_boundaries!(state, variables, turbulencesetup)

    return
end

"""
    set_vertical_boundaries!(state, variables::BoundaryGWIntegrals, boundaries)

Set vertical boundaries for GW integral fields. Dispatches based on WKB mode.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWIntegrals,
    boundaries::AbstractBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    (; tracersetup) = state.namelists.tracer 

    set_vertical_boundaries!(state, variables, boundaries, wkb_mode)
    set_tracer_vertical_boundaries!(state, variables, boundaries, wkb_mode, tracersetup)
    return
end

"""
    set_vertical_boundaries!(state, variables::BoundaryGWIntegrals, boundaries, wkb_mode::AbstractWKBMode)

Set vertical boundaries for basic GW integral fields (uw, vw, e) with symmetric conditions.
"""
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
    set_vertical_boundaries!(state, variables::BoundaryGWIntegrals, boundaries, wkb_mode::MultiColumn)

Set vertical boundaries for extended GW integral fields in multi-column mode.
"""
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
    set_vertical_boundaries!(state, variables::BoundaryGWTendencies, boundaries)

Set vertical boundaries for GW tendency fields. Dispatches based on WKB mode.
"""
function set_vertical_boundaries!(
    state::State,
    variables::BoundaryGWTendencies,
    boundaries::AbstractBoundaries,
)
    (; wkb_mode) = state.namelists.wkb
    (; tracersetup) = state.namelists.tracer
    set_vertical_boundaries!(state, variables, boundaries, wkb_mode)
    set_tracer_vertical_boundaries!(state, variables, boundaries, wkb_mode, tracersetup)
    return
end

"""
    set_vertical_boundaries!(state, variables::BoundaryGWTendencies, boundaries, wkb_mode::AbstractWKBMode)

Set vertical boundaries for basic GW tendency fields (dudt, dvdt) with symmetric conditions.
"""
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

"""
    set_vertical_boundaries!(state, variables::BoundaryGWTendencies, boundaries, wkb_mode::MultiColumn)

Set vertical boundaries for GW tendency fields in multi-column mode, including dthetadt.
"""
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

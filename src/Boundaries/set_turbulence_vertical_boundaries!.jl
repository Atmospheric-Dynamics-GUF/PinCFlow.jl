"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::NoTurbulence,
)
    return
end

"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; turbulencepredictands) = state.turbulence

    for field in fieldnames(TurbulencePredictands)
        set_vertical_boundaries_of_field!(
            getfield(turbulencepredictands, field),
            namelists,
            domain,
            zboundaries,
            -,
        )
    end

    return
end

"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::NoTurbulence,
)
    return
end

"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryReconstructions,
    turbulencesetup::AbstractTurbulence,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; turbulencereconstructions) = state.turbulence

    for field in fieldnames(TurbulenceReconstructions)
        set_vertical_boundaries_of_field!(
            getfield(turbulencereconstructions, field),
            namelists,
            domain,
            zboundaries,
        )
    end

    return
end

"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::NoTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::NoTurbulence,
)
    return
end

"""
```julia
set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::AbstractTurbulence,
)
```
"""
function set_turbulence_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    turbulencesetup::AbstractTurbulence,
)
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; turbulencefluxes) = state.turbulence

    if ko == 0
        for field in fieldnames(TurbulenceFluxes)
            getfield(turbulencefluxes, field)[:, :, k0 - 1, 3] .= 0.0
        end
    end

    if ko + nzz == sizezz
        for field in fieldnames(TurbulenceFluxes)
            getfield(turbulencefluxes, field)[:, :, k1, 3] .= 0.0
        end
    end

    return
end

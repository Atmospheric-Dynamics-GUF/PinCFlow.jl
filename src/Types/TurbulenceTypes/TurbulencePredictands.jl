"""
```julia
TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    tke::A
    tte::A
end

"""
```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
```
"""
function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    (; turbulencesetup) = namelists.turbulence

    return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        turbulencesetup,
        variables,
    )
end

"""
```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::NoTurbulence,
    variables::Variables,
)
```
"""
function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::NoTurbulence,
    variables::Variables,
)
    tke = zeros(0, 0, 0)
    tte = zeros(0, 0, 0)

    return TurbulencePredictands(tke, tte)
end

"""
```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::AbstractTurbulence,
    variables::Variables,
)
```
"""
function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::AbstractTurbulence,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands
    (; lref, tref) = constants

    tke =
        ones(nxx, nyy, nzz) .* 0.1 .* (tref ^ 2.0) ./ (lref ^ 2.0) .*
        (rho .+ rhostrattfc)
    tte =
        ones(nxx, nyy, nzz) .* 0.1 .* (tref ^ 2.0) ./ (lref ^ 2.0) .*
        (rho .+ rhostrattfc)

    return TurbulencePredictands(tke, tte)
end

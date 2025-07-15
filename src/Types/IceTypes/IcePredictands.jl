"""
```julia
IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```
"""
struct IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    n::A # ice crystal number concentration n
    q::A # ice mixing ratio q
    qv::A # vapor mixing ratio qv
end

"""
```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
```
"""
function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    (; icesetup) = namelists.ice

    return IcePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        icesetup,
        variables,
    )
end

"""
```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::NoIce,
    variables::Variables,
)
```
"""
function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::NoIce,
    variables::Variables,
)
    n = zeros(0, 0, 0)
    q = zeros(0, 0, 0)
    qv = zeros(0, 0, 0)

    return IcePredictands(n, q, qv)
end

"""
```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::AbstractIce,
    variables::Variables,
)
```
"""
function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::AbstractIce,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands

    n = zeros(nxx, nyy, nzz)
    q = zeros(nxx, nyy, nzz)
    qv = zeros(nxx, nyy, nzz)

    n .= ztfc .* (rho .+ rhostrattfc)
    q .= ztfc .* (rho .+ rhostrattfc)
    qv .= ztfc .* (rho .+ rhostrattfc)

    return IcePredictands(n, q, qv)
end

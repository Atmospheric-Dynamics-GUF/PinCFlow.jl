struct TracerPredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    chi::A
    # can be extended to included more passive tracers
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    (; tracersetup) = namelists.tracer

    return TracerPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        tracersetup,
        variables,
    )
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracersetup::AbstractTracer,
    variables::Variables,
)
    chi = zeros(0, 0, 0)

    return TracerPredictands(chi)
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracersetup::LinearTracer,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands
    (; testcase) = namelists.setting
    (; lref) = constants

    alphatracer = lref
    chi = zeros(nxx, nyy, nzz)
    chi .= alphatracer .* ztfc

    initialize_tracer_wave_packet!(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
        alphatracer,
        chi,
        testcase,
    )

    chi .= chi .* (rho .+ rhostrattfc)

    return TracerPredictands(chi)
end

function TracerPredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    tracersetup::LikeDensity,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands
    (; testcase) = namelists.setting
    (; lref) = constants

    chi = zeros(nxx, nyy, nzz)

    chi .= rho .+ rhostrattfc

    return TracerPredictands(chi)
end

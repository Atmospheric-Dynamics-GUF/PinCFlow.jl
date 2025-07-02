struct Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceTendencies,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
}
    turbulencepredictands::A
    turbulencetendencies::B
    turbulenceauxiliaries::C
    turbulencereconstructions::D
    turbulencefluxes::E
end

function Turbulence(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    turbulencepredictands = TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    turbulencetendencies = TurbulenceTendencies(namelists, domain)
    turbulenceauxiliaries = TurbulenceAuxiliaries(constants)
    turbulencereconstructions = TurbulenceReconstructions(namelists, domain)
    turbulencefluxes = TurbulenceFluxes(namelists, domain)

    return Turbulence(
        turbulencepredictands,
        turbulencetendencies,
        turbulenceauxiliaries,
        turbulencereconstructions,
        turbulencefluxes,
    )
end

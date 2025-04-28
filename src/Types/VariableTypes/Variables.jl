struct Variables{
    A <: Predictands,
    B <: Tendencies,
    C <: Backups,
    D <: Auxiliaries,
    E <: Reconstructions,
    F <: Fluxes,
}
    predictands::A
    tendencies::B
    backups::C
    auxiliaries::D
    reconstructions::E
    fluxes::F
end

function Variables(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)

    # Initialize all fields.
    predictands = Predictands(namelists, constants, domain, atmosphere)
    tendencies = Tendencies(namelists, domain)
    backups = Backups(domain)
    auxiliaries = Auxiliaries(domain)
    reconstructions = Reconstructions(domain)
    fluxes = Fluxes(namelists, domain)

    # Return a Variables instance.
    return Variables(
        predictands,
        tendencies,
        backups,
        auxiliaries,
        reconstructions,
        fluxes,
    )
end

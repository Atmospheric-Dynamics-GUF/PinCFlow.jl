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

function Variables(namelists::Namelists, constants::Constants, domain::Domain)

  # Initialize all fields.
  predictands = Predictands(namelists, constants, domain)
  tendencies = Tendencies(namelists, domain)
  backups = Backups(namelists, domain)
  auxiliaries = Auxiliaries(namelists, domain)
  reconstructions = Reconstructions(namelists, domain)
  fluxes = Fluxes(domain)

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

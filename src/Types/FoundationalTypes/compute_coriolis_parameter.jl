function compute_coriolis_parameter(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    coriolis_mode::FPlane,
)
    (; coriolis_frequency) = namelists.atmosphere
    (; tref) = constants
    (; nyy) = domain
    return coriolis_frequency .* tref .* ones(nyy)
end

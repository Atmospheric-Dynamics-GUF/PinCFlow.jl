struct Predictands{A<:OffsetArray{<:AbstractFloat,3}}
    rho::A
    rhop::A
    u::A
    v::A
    w::A
    pip::A
end

struct Tendencies{A<:OffsetArray{<:AbstractFloat,3}}
    drho::A
    drhop::A
    du::A
    dv::A
    dw::A
    dpip::A
end

struct Backups{A<:OffsetArray{<:AbstractFloat,3}}
    rhoold::A
    rhopold::A
    uold::A
    vold::A
    wold::A
end

struct Auxiliaries{A<:OffsetArray{<:AbstractFloat,3}}
    rhobar::A
    rhopbar::A
    ubar::A
    vbar::A
    wbar::A
end

struct Reconstructions{A<:OffsetArray{<:AbstractFloat,5}}
    rhotilde::A
    rhoptilde::A
    utilde::A
    vtilde::A
    wtilde::A
end

struct Fluxes{A<:OffsetArray{<:AbstractFloat,4}}
    phirho::A
    phirhop::A
    phiu::A
    phiv::A
    phiw::A
end

struct Variables{
    A<:Predictands,
    B<:Tendencies,
    C<:Backups,
    D<:Auxiliaries,
    E<:Reconstructions,
    F<:Fluxes,
}
    predictands::A
    tendencies::B
    backups::C
    auxiliaries::D
    reconstructions::E
    fluxes::F
end

function Predictands(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  model::PseudoIncompressible,
  testcase::MountainWave,
)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; uref) = constants
  (; nx, ny, nz) = domain

  # Initialize the predictands.
  (rho, rhop, u, v, w, pip) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:6
  )

  # Initial winds.
  u .= backgroundflow_dim[1] / uref
  v .= backgroundflow_dim[2] / uref
  w .= backgroundflow_dim[3] / uref

  # Return a Predictands instance.
  return Predictands(rho, rhop, u, v, w, pip)
end

function Tendencies(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the tendencies.
  (drho, drhop, du, dv, dw, dpip) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:6
  )

  # Return a Variables instance.
  return Tendencies(drho, drhop, du, dv, dw, dpip)
end

function Backups(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the backups.
  (rhoold, rhopold, uold, vold, wold) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:5
  )

  # Return a Backups instance.
  return Backups(rhoold, rhopold, uold, vold, wold)
end

function Auxiliaries(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the auxiliaries.
  (rhobar, rhopbar, ubar, vbar, wbar) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
    ) for i in 1:5
  )

  # Return an Auxiliaries instance.
  return Auxiliaries(rhobar, rhopbar, ubar, vbar, wbar)
end

function Reconstructions(namelists::Namelists, domain::Domain)

  # Get parameters.
  (; nbx, nby, nbz) = namelists.domain
  (; nx, ny, nz) = domain

  # Initialize the reconstructed variables.
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) = (
    OffsetArray(
      zeros((nx + 2 * nbx + 1, ny + 2 * nby + 1, nz + 2 * nbz + 1, 3, 2)),
      (-nbx):(nx + nbx),
      (-nby):(ny + nby),
      (-nbz):(nz + nbz),
      1:3,
      0:1,
    ) for i in 1:5
  )

  # Return a Reconstructions instance.
  return Reconstructions(rhotilde, rhoptilde, utilde, vtilde, wtilde)
end

function Fluxes(domain::Domain)

  # Get parameters.
  (; nx, ny, nz) = domain

  # Initialize the fluxes.
  (phirho, phirhop, phiu, phiv, phiw) = (
    OffsetArray(
      zeros((nx + 2, ny + 2, nz + 2, 3)),
      (-1):(nx),
      (-1):(ny),
      (-1):(nz),
      1:3,
    ) for i in 1:5
  )

  # Return a Fluxes instance.
  return Fluxes(phirho, phirhop, phiu, phiv, phiw)
end

function Variables(namelists::Namelists, constants::Constants, domain::Domain)

  # Get parameters.
  (; model, testcase) = namelists.setting

  # Initialize all fields.
  predictands = Predictands(namelists, constants, domain, model, testcase)
  tendencies = Tendencies(namelists, domain)
  backups = Backups(namelists, domain)
  auxiliaries = Auxiliaries(namelists, domain)
  reconstructions = Reconstructions(namelists, domain)
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

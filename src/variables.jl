using OffsetArrays
"""
Atmospheric variables.
"""
struct Variables{
    B<:Bool,
    A3OF<:OffsetArray{<:AbstractFloat,3},
}

    # Bools for prognostic equations and Poisson problem.
    updatemass::B
    predictmomentum::B
    correctmomentum::B

    # Prognostic variables.
    # as named tuple and with symbols :rho, :rhop, :u, :v, :w, :pip and fields
    prognostic_fields::NamedTuple
    # Tendencies.
    # same as for prognostic variables
    tendencies::NamedTuple
    # # Pressure correction.
    dpip::A3OF
    # # Saved variables.
    history::NamedTuple
    # # Reconstructed variables.
    reconstructed::NamedTuple
end

function Variables(pars::Parameters)
    prognostic = prognostic_fields(pars.domain)
    tendency = tendency_fields(pars.domain)
    reconstruced = reconstructed_fields(pars.domain)
    history = history_fields(prognostic)
    dpip = copy(prognostic.u)
    # TODO: default values for boolean flags?
    Variables(false, false, false, prognostic, tendency, dpip, history, reconstruced)
end

function history_fields(prognostic_fields)
    return NamedTuple{keys(prognostic_fields)}((copy(f) for f in prognostic_fields))
end

function reconstructed_fields(domain::DomainParameters)
    nx, ny, nz = domain.sizex, domain.sizey, domain.sizez
    nbx, nby, nbz = domain.nbx, domain.nby, domain.nbz
    ndim = 3
    # TODO: factor this in own function
    f = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz, ndim,
            2),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz),
        1:ndim,
        0:1)
    return NamedTuple{(:rho, :rhop, :u, :v, :w)}((f,
        copy(f), copy(f), copy(f), copy(f)))
end

function prognostic_fields(domain::DomainParameters)
    nx, ny, nz = domain.sizex, domain.sizey, domain.sizez
    nbx, nby, nbz = domain.nbx, domain.nby, domain.nbz
    u = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz))
    # TODO: pip == exner?
    v, w, rho, pip, rhop = (copy(u), copy(u), copy(u), copy(u), copy(u))

    return NamedTuple{(:u, :v, :w, :rho, :pip, :rhop)}((u, v, w, rho, pip, rhop))
end

function tendency_fields(domain::DomainParameters)
    nx, ny, nz = domain.sizex, domain.sizey, domain.sizez
    nbx, nby, nbz = domain.nbx, domain.nby, domain.nbz
    ndim = 3
    drho = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz))
    drhop = copy(drho)
    dmom = OffsetArray(zeros(Float64, nx + 1 + 2nbx, ny + 1 + 2nby, nz + 1 + 2nbz, ndim),
        (-nbx):(nx+nbx),
        (-nby):(ny+nby),
        (-nbz):(nz+nbz),
        1:ndim)
    return NamedTuple{(:drho, :drhop, :dmom)}((drho, drhop, dmom))
end

"""
```julia
DomainNamelist{
    A <: Integer,
    B <: Integer,
    C <: Integer,
    D <: Integer,
    E <: Integer,
    F <: Integer,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: Integer,
    K <: Integer,
    L <: Integer,
    M <: MPI.Comm,
}
```

Namelist for parameters describing the model domain.

```julia
DomainNamelist(;
    x_size::Integer = 3,
    y_size::Integer = 3,
    z_size::Integer = 3,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx::AbstractFloat = 1.0E+3,
    ly::AbstractFloat = 1.0E+3,
    lz::AbstractFloat = 1.0E+3,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
```

Construct a `DomainNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `x_size::A`: Number of grid cells in ``\\widehat{x}``-direction.

  - `y_size::B`: Number of grid cells in ``\\widehat{y}``-direction.

  - `z_size::C`: Number of grid cells in ``\\widehat{z}``-direction.

  - `nbx::D`: Number of boundary/halo cells in ``\\widehat{x}``-direction.

  - `nby::E`: Number of boundary/halo cells in ``\\widehat{y}``-direction.

  - `nbz::F`: Number of boundary/halo cells in ``\\widehat{z}``-direction.

  - `lx::G`: Domain extent in ``\\widehat{x}``-direction.

  - `ly::H`: Domain extent in ``\\widehat{y}``-direction.

  - `lz::I`: Domain extent in ``\\widehat{z}``-direction.

  - `npx::J`: Number of MPI processes in ``\\widehat{x}``-direction.

  - `npy::K`: Number of MPI processes in ``\\widehat{y}``-direction.

  - `npz::L`: Number of MPI processes in ``\\widehat{z}``-direction.

  - `base_comm::M`: MPI base communicator.
"""
struct DomainNamelist{
    A <: Integer,
    B <: Integer,
    C <: Integer,
    D <: Integer,
    E <: Integer,
    F <: Integer,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: Integer,
    K <: Integer,
    L <: Integer,
    M <: MPI.Comm,
}
    x_size::A
    y_size::B
    z_size::C
    nbx::D
    nby::E
    nbz::F
    lx::G
    ly::H
    lz::I
    npx::J
    npy::K
    npz::L
    base_comm::M
end

function DomainNamelist(;
    x_size::Integer = 3,
    y_size::Integer = 3,
    z_size::Integer = 3,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx::AbstractFloat = 1.0E+3,
    ly::AbstractFloat = 1.0E+3,
    lz::AbstractFloat = 1.0E+3,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
    return DomainNamelist(
        x_size,
        y_size,
        z_size,
        nbx,
        nby,
        nbz,
        lx,
        ly,
        lz,
        npx,
        npy,
        npz,
        base_comm,
    )
end

"""
```julia
DomainNamelist{A <: Integer, B <: NTuple{2, <:AbstractFloat}, C <: MPI.Comm}
```

Namelist for parameters describing the model domain.

```julia
DomainNamelist(;
    ndx::Integer = 3,
    ndy::Integer = 3,
    ndz::Integer = 3,
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

  - `ndx::A`: Number of grid cells in ``\\widehat{x}``-direction.

  - `ndy::A`: Number of grid cells in ``\\widehat{y}``-direction.

  - `ndz::A`: Number of grid cells in ``\\widehat{z}``-direction.

  - `nbx::A`: Number of boundary/halo cells in ``\\widehat{x}``-direction.

  - `nby::A`: Number of boundary/halo cells in ``\\widehat{y}``-direction.

  - `nbz::A`: Number of boundary/halo cells in ``\\widehat{z}``-direction.

  - `lx::B`: Domain extent in ``\\widehat{x}``-direction.

  - `ly::B`: Domain extent in ``\\widehat{y}``-direction.

  - `lz::B`: Domain extent in ``\\widehat{z}``-direction.

  - `npx::A`: Number of MPI processes in ``\\widehat{x}``-direction.

  - `npy::A`: Number of MPI processes in ``\\widehat{y}``-direction.

  - `npz::A`: Number of MPI processes in ``\\widehat{z}``-direction.

  - `base_comm::C`: MPI base communicator.
"""
struct DomainNamelist{A <: Integer, B <: AbstractFloat, C <: MPI.Comm}
    ndx::A
    ndy::A
    ndz::A
    nbx::A
    nby::A
    nbz::A
    lx::B
    ly::B
    lz::B
    npx::A
    npy::A
    npz::A
    base_comm::C
end

function DomainNamelist(;
    ndx::Integer = 3,
    ndy::Integer = 3,
    ndz::Integer = 3,
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
        ndx,
        ndy,
        ndz,
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

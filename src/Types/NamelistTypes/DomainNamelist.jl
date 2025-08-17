"""
```julia
DomainNamelist{A <: Integer, B <: NTuple{2, <:AbstractFloat}, C <: MPI.Comm}
```

Namelist for parameters describing the model domain.

```julia
DomainNamelist(;
    sizex::Integer = 3,
    sizey::Integer = 3,
    sizez::Integer = 3,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    ly_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    lz_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
```

Construct a `DomainNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `sizex::A`: Number of grid cells in ``\\widehat{x}``-direction.

  - `sizey::A`: Number of grid cells in ``\\widehat{y}``-direction.

  - `sizez::A`: Number of grid cells in ``\\widehat{z}``-direction.

  - `nbx::A`: Number of boundary/halo cells in ``\\widehat{x}``-direction.

  - `nby::A`: Number of boundary/halo cells in ``\\widehat{y}``-direction.

  - `nbz::A`: Number of boundary/halo cells in ``\\widehat{z}``-direction.

  - `lx_dim::B`: Domain boundaries in ``\\widehat{x}``-direction.

  - `ly_dim::B`: Domain boundaries in ``\\widehat{y}``-direction.

  - `lz_dim::B`: Domain boundaries in ``\\widehat{z}``-direction.

  - `npx::A`: Number of MPI processes in ``\\widehat{x}``-direction.

  - `npy::A`: Number of MPI processes in ``\\widehat{y}``-direction.

  - `npz::A`: Number of MPI processes in ``\\widehat{z}``-direction.

  - `base_comm::C`: MPI base communicator.
"""
struct DomainNamelist{
    A <: Integer,
    B <: NTuple{2, <:AbstractFloat},
    C <: MPI.Comm,
}
    sizex::A
    sizey::A
    sizez::A
    nbx::A
    nby::A
    nbz::A
    lx_dim::B
    ly_dim::B
    lz_dim::B
    npx::A
    npy::A
    npz::A
    base_comm::C
end

function DomainNamelist(;
    sizex::Integer = 3,
    sizey::Integer = 3,
    sizez::Integer = 3,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    ly_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    lz_dim::NTuple{2, <:AbstractFloat} = (0.0E+0, 1.0E+3),
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
    return DomainNamelist(
        sizex,
        sizey,
        sizez,
        nbx,
        nby,
        nbz,
        lx_dim,
        ly_dim,
        lz_dim,
        npx,
        npy,
        npz,
        base_comm,
    )
end

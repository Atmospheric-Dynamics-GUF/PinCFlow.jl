"""
```julia
DomainNamelist{A <: Int, B <: Float64, C <: MPI.Comm}
```

Namelist for parameters describing the model domain.

```julia
DomainNamelist(;
    x_size::Integer = 1,
    y_size::Integer = 1,
    z_size::Integer = 1,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx::Real = 1.0E+3,
    ly::Real = 1.0E+3,
    lz::Real = 1.0E+3,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
```

Construct a `DomainNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `x_size::A`: Number of grid cells in ``\\hat{x}``-direction.

  - `y_size::A`: Number of grid cells in ``\\hat{y}``-direction.

  - `z_size::A`: Number of grid cells in ``\\hat{z}``-direction.

  - `nbx::A`: Number of boundary/halo cells in ``\\hat{x}``-direction.

  - `nby::A`: Number of boundary/halo cells in ``\\hat{y}``-direction.

  - `nbz::A`: Number of boundary/halo cells in ``\\hat{z}``-direction.

  - `lx::B`: Domain extent in ``\\hat{x}``-direction.

  - `ly::B`: Domain extent in ``\\hat{y}``-direction.

  - `lz::B`: Domain extent in ``\\hat{z}``-direction.

  - `npx::A`: Number of MPI processes in ``\\hat{x}``-direction.

  - `npy::A`: Number of MPI processes in ``\\hat{y}``-direction.

  - `npz::A`: Number of MPI processes in ``\\hat{z}``-direction.

  - `base_comm::C`: MPI base communicator.
"""
struct DomainNamelist{A <: Int, B <: Float64, C <: MPI.Comm}
    x_size::A
    y_size::A
    z_size::A
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
    x_size::Integer = 1,
    y_size::Integer = 1,
    z_size::Integer = 1,
    nbx::Integer = 3,
    nby::Integer = 3,
    nbz::Integer = 3,
    lx::Real = 1.0E+3,
    ly::Real = 1.0E+3,
    lz::Real = 1.0E+3,
    npx::Integer = 1,
    npy::Integer = 1,
    npz::Integer = 1,
    base_comm::MPI.Comm = MPI.COMM_WORLD,
)::DomainNamelist
    return DomainNamelist(
        Int(x_size),
        Int(y_size),
        Int(z_size),
        Int(nbx),
        Int(nby),
        Int(nbz),
        Float64(lx),
        Float64(ly),
        Float64(lz),
        Int(npx),
        Int(npy),
        Int(npz),
        base_comm,
    )
end

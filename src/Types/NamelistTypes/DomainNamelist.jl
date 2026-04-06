"""
```julia
DomainNamelist
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

  - `x_size::Int`: Number of grid cells in ``\\hat{x}``-direction.

  - `y_size::Int`: Number of grid cells in ``\\hat{y}``-direction.

  - `z_size::Int`: Number of grid cells in ``\\hat{z}``-direction.

  - `nbx::Int`: Number of boundary/halo cells in ``\\hat{x}``-direction.

  - `nby::Int`: Number of boundary/halo cells in ``\\hat{y}``-direction.

  - `nbz::Int`: Number of boundary/halo cells in ``\\hat{z}``-direction.

  - `lx::Float64`: Domain extent in ``\\hat{x}``-direction.

  - `ly::Float64`: Domain extent in ``\\hat{y}``-direction.

  - `lz::Float64`: Domain extent in ``\\hat{z}``-direction.

  - `npx::Int`: Number of MPI processes in ``\\hat{x}``-direction.

  - `npy::Int`: Number of MPI processes in ``\\hat{y}``-direction.

  - `npz::Int`: Number of MPI processes in ``\\hat{z}``-direction.

  - `base_comm::C`: MPI base communicator.
"""
struct DomainNamelist
    x_size::Int
    y_size::Int
    z_size::Int
    nbx::Int
    nby::Int
    nbz::Int
    lx::Float64
    ly::Float64
    lz::Float64
    npx::Int
    npy::Int
    npz::Int
    base_comm::MPI.Comm
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

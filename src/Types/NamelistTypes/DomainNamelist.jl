"""
```julia
DomainNamelist{A <: Integer, B <: NTuple{2, <:AbstractFloat}, C <: MPI.Comm}
```

Namelist for the model domain (see constructor for parameter descriptions).
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

"""
```julia
DomainNamelist(;
    sizex = 4,
    sizey = 4,
    sizez = 4,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0E+0, 1.0E+3),
    ly_dim = (0.0E+0, 1.0E+3),
    lz_dim = (0.0E+0, 1.0E+3),
    npx = 1,
    npy = 1,
    npz = 1,
    base_comm = MPI.COMM_WORLD,
)
```

Construct a DomainNamelist, which holds parameters for the spatial domain.

# Arguments

  - `sizex`: Number of grid cells in the x-direction.
  - `sizey`: Number of grid cells in the y-direction.
  - `sizez`: Number of grid cells in the z-direction.
  - `nbx`: Number of boundary cells in the x-direction.
  - `nby`: Number of boundary cells in the y-direction.
  - `nbz`: Number of boundary cells in the z-direction.
  - `lx_dim`: Dimensions of the domain in the x-direction.
  - `ly_dim`: Dimensions of the domain in the y-direction.
  - `lz_dim`: Dimensions of the domain in the z-direction.
  - `npx`: Number of processors in the x-direction.
  - `npy`: Number of processors in the y-direction.
  - `npz`: Number of processors in the z-direction.
  - `base_comm`: MPI base communicator.
"""
function DomainNamelist(;
    sizex = 4,
    sizey = 4,
    sizez = 4,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0E+0, 1.0E+3),
    ly_dim = (0.0E+0, 1.0E+3),
    lz_dim = (0.0E+0, 1.0E+3),
    npx = 1,
    npy = 1,
    npz = 1,
    base_comm = MPI.COMM_WORLD,
)
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

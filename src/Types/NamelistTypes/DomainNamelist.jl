struct DomainNamelist{A <: Integer, B <: NTuple{2, <:AbstractFloat}}
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
end

"""
DomainNamelist(; <keyword arguments>)

Construct a DomainNamelist, which holds parameters for the spatial domain.

# Arguments

  - `sizex::Integer = 4`: Number of grid cells in the x-direction.
  - `sizey::Integer = 4`: Number of grid cells in the y-direction.
  - `sizez::Integer = 4`: Number of grid cells in the z-direction.
  - `nbx::Integer = 3`: Number of boundary cells in the x-direction.
  - `nby::Integer = 3`: Number of boundary cells in the y-direction.
  - `nbz::Integer = 3`: Number of boundary cells in the z-direction.
  - `lx_dim::NTuple{2, <:AbstractFloat} = (0.0, 1.0E+3)`: Dimensions of the domain in the x-direction.
  - `ly_dim::NTuple{2, <:AbstractFloat} = (0.0, 1.0E+3)`: Dimensions of the domain in the y-direction.
  - `lz_dim::NTuple{2, <:AbstractFloat} = (0.0, 1.0E+3)`: Dimensions of the domain in the z-direction.
  - `npx::Integer = 1`: Number of processors in the x-direction.
  - `npy::Integer = 1`: Number of processors in the y-direction.
  - `npz::Integer = 1`: Number of processors in the z-direction.
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
    )
end

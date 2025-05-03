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

function DomainNamelist(;
    sizex = 4,
    sizey = 4,
    sizez = 4,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 1.0E+3),
    ly_dim = (0.0, 1.0E+3),
    lz_dim = (0.0, 1.0E+3),
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

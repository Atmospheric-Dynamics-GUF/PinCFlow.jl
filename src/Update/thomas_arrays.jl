function thomas_arrays end

function thomas_arrays(state::State, nx::Integer, ny::Integer, nz::Integer)

    (; athglob, bthglob, cthglob, fthglob, qthglob, pthglob, qthglob_bc, fthglob_bc) =
        state.variables.auxiliaries
    
    @ivy ath = athglob[1:nx, 1:ny, 1:nz]
    @ivy bth = bthglob[1:nx, 1:ny, 1:nz]
    @ivy cth = cthglob[1:nx, 1:ny, 1:nz]
    @ivy fth = fthglob[1:nx, 1:ny, 1:nz]
    @ivy qth = qthglob[1:nx, 1:ny, 1:nz]
    @ivy pth = pthglob[1:nx, 1:ny]
    @ivy fth_bc = fthglob_bc[1:nx, 1:ny]
    @ivy qth_bc = qthglob_bc[1:nx, 1:ny]

    return ath, bth, cth, fth, qth, pth, fth_bc, qth_bc 
end
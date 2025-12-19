struct SubGrid{A <: AbstractFloat, B <: Integer}

	sizex2::B
	sizey2::B
	sizez2::B

	dxsc::A
	dysc::A
	dzsc::A

	nxnscx::B
	nynscy::B
	nznscz::B

	nxnscxx::B
	nynscyy::B
	nznsczz::B

	i02::B
	j02::B
	k02::B
	i12::B
	j12::B
	k12::B

    x2::AbstractVector{A}
    y2::AbstractVector{A}
    z2tfc::AbstractArray{A, 3}

end

function SubGrid(namelists::Namelists, domain::Domain, grid::Grid)

	(; cloudcover) = namelists.ice
	return SubGrid(namelists, domain, grid, cloudcover)
end

function SubGrid(namelists::Namelists,
	domain::Domain,
	grid::Grid,
	cloudcover::CloudCoverOff,
)

	(; nbx, nby, nbz) = namelists.domain
	(; x_size, y_size, z_size) = namelists.domain
	(; nx, ny, nz, i0, j0, k0, i1, j1, k1) = domain
	(; dx, dy, dz) = grid

	sizex2 = x_size
	sizey2 = y_size
	sizez2 = z_size

	dxsc = dx
	dysc = dy
	dzsc = dz

	nxnscx = nx
	nynscy = ny
	nznscz = nz

	nxnscxx = nxnscx + 2 * nbx
	nynscyy = nynscy + 2 * nby
	nznsczz = nznscz + 2 * nbz

	i02 = i0
	j02 = j0
	k02 = k0
	i12 = i1
	j12 = j1
	k12 = k1

	x2 = zeros(nxnscxx)
	y2 = zeros(nynscyy)
	z2tfc = zeros(nxnscxx, nynscyy, nznsczz)

	return SubGrid(sizex2, sizey2, sizez2,
		dxsc, dysc, dzsc,
		nxnscx, nynscy, nznscz,
		nxnscxx, nynscyy, nznsczz, 
		i02, j02, k02, 
		i12, j12, k12, 
		x2, y2, z2tfc)

end

function SubGrid(namelists::Namelists,
	domain::Domain,
	grid::Grid,
	cloudcover::CloudCoverOn,
)

	(; nscx, nscy, nscz) = namelists.ice
	(; nbx, nby, nbz) = namelists.domain
	(; x_size, y_size, z_size) = namelists.domain
	(; lx, ly, dx, dy, dz) = grid
    (; x, y, zctilde, jac) = grid
	(; nx, ny, nz, io, jo, ko,
	i0, j0, k0, i1, j1, k1) = domain

		dxsc = dx / nscx
		dysc = dy / nscy
		dzsc = dz / nscz

		nxnscx = nx * nscx
		nynscy = ny * nscy
		nznscz = nz * nscz

		nxnscxx = nxnscx + 2 * nbx * nscx
		nynscyy = nynscy + 2 * nby * nscy
		nznsczz = nznscz + 2 * nbz * nscz

		i02 = nbx * nscx + 1
		j02 = nby * nscy + 1
		k02 = nbz * nscz + 1	

		i12 = i02 + nxnscx - 1
		j12 = j02 + nynscy - 1
		k12 = k02 + nznscz - 1

		sizex2 = x_size * nscx
		sizey2 = y_size * nscy
		sizez2 = z_size * nscz

		# Compute sub grid x-coordinate.

		x2 = zeros(nxnscxx)
		for ix in 1:nxnscxx
			x2[ix] = -lx/2. + (ix + io*nscx - i02) * dxsc + dxsc / 2
		end

		# Compute y-coordinate.

		y2 = zeros(nynscyy)
		for jy in 1:nynscyy
			y2[jy] = -ly/2. + (jy + jo*nscy - j02) * dysc + dysc / 2
		end

		# Initialize the physical layers.
		(z2tfc) = zeros(nxnscxx, nynscyy, nznsczz)

		for ix in i0:i1
			for ii in 1:nscx
				ii2 = (ix - 1) * nscx + ii
				for jy in j0:j1
					for jj in 1:nscy
						jj2 = (jy - 1) * nscy + jj
						for kz in k0:k1
							for kk in 1:nscz
								kk2 = (kz - 1) * nscz + kk

								#subcell center
								#xsc = x[ix+io] - dx / 2.0 +
								#	  (ii - 0.5) * dxsc
								#ysc = y[jy+jo] - dy / 2.0 +
								#	  (jj - 0.5) * dysc
								zsc = zctilde[ix, jy, kz-1] +
									  (kk - 0.5) * dzsc * jac[ix, jy, kz]

								z2tfc[ii2, jj2, kk2] = zsc

							end
						end
					end
				end
			end
		end

	return SubGrid(sizex2, sizey2, sizez2,
		dxsc, dysc, dzsc,
		nxnscx, nynscy, nznscz,
		nxnscxx, nynscyy, nznsczz, 
		i02, j02, k02, 
		i12, j12, k12, 
		x2, y2, z2tfc)

end

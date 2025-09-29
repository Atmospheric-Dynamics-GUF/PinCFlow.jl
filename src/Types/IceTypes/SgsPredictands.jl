struct SgsPredictands{A <: AbstractArray{<:AbstractFloat, 3}}
	n::A # ice crystal number concentration n
	q::A # ice mixing ratio q
	qv::A # vapor mixing ratio qv
end

function SgsPredictands(
	namelists::Namelists,
	domain::Domain,
	icepredictands::IcePredictands,
	subgrid::SubGrid,
)
	(; cloudcover) = namelists.ice

	return SgsPredictands(
		namelists,
		domain,
		icepredictands,
		subgrid,
		cloudcover,
	)
end

function SgsPredictands(
	namelists::Namelists,
	domain::Domain,
	icepredictands::IcePredictands,
	subgrid::SubGrid,
	cloudcover::CloudCoverOff,
)
	n = zeros(0, 0, 0)
	q = zeros(0, 0, 0)
	qv = zeros(0, 0, 0)

	return SgsPredictands(n, q, qv)
end

function SgsPredictands(
	namelists::Namelists,
	domain::Domain,
	icepredictands::IcePredictands,
	subgrid::SubGrid,
	cloudcover::CloudCoverOn,
)

	(; i0, i1, j0, j1, k0, k1) = domain
	(; nscx, nscy, nscz) = namelists.ice
	(; n, q, qv) = icepredictands
	(; nxnscxx, nynscyy, nznsczz) = subgrid

	sgs_n = zeros(nxnscxx, nynscyy, nznsczz)
	sgs_q = zeros(nxnscxx, nynscyy, nznsczz)
	sgs_qv = zeros(nxnscxx, nynscyy, nznsczz)

	# initialize
	for ix in i0:i1, jy in j0:j1, kz in k0:k1
		sgs_n[((ix-1)*nscx+1):(ix*nscx), ((jy-1)*nscy+1):(jy*nscy), ((kz-1)*nscz+1):(kz*nscz)] .= n[ix, jy, kz]
		sgs_q[((ix-1)*nscx+1):(ix*nscx), ((jy-1)*nscy+1):(jy*nscy), ((kz-1)*nscz+1):(kz*nscz)] .= q[ix, jy, kz]
		sgs_qv[((ix-1)*nscx+1):(ix*nscx), ((jy-1)*nscy+1):(jy*nscy), ((kz-1)*nscz+1):(kz*nscz)] .= qv[ix, jy, kz]
	end

	return SgsPredictands(sgs_n, sgs_q, sgs_qv)

end

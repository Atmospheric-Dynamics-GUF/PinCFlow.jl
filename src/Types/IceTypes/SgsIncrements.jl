struct SgsIncrements{A <: AbstractArray{<:AbstractFloat, 3}}
	n::A # ice crystal number concentration n
	q::A # ice mixing ratio q
	qv::A # vapor mixing ratio qv
end

function SgsIncrements(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
)
	(; icesetup) = namelists.ice

	return SgsIncrements(
		namelists,
		domain,
		subgrid,
		icesetup,
	)
end

function SgsIncrements(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
	icesetup::NoIce,
)
	n = zeros(0, 0, 0)
	q = zeros(0, 0, 0)
	qv = zeros(0, 0, 0)

	return SgsIncrements(n, q, qv)
end

function SgsIncrements(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
)
	(; cloudcover) = namelists.ice

	SgsIncrements(
		namelists,
		domain,
		subgrid,
		cloudcover,
	)
end

function SgsIncrements(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
	cloudcover::CloudCoverOff,
)

	sgs_n = zeros(0, 0, 0)
	sgs_q = zeros(0, 0, 0)
	sgs_qv = zeros(0, 0, 0)

	return SgsIncrements(sgs_n, sgs_q, sgs_qv)
end

function SgsIncrements(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
	cloudcover::CloudCoverOn,
)
	#(; i0, i1, j0, j1, k0, k1) = domain
	#(; nscx, nscy, nscz) = namelists.ice
	#(; n, q, qv) = icepredictands
	(; nxnscxx, nynscyy, nznsczz) = subgrid

	sgs_n = zeros(nxnscxx, nynscyy, nznsczz)
	sgs_q = zeros(nxnscxx, nynscyy, nznsczz)
	sgs_qv = zeros(nxnscxx, nynscyy, nznsczz)

	return SgsIncrements(sgs_n, sgs_q, sgs_qv)
end

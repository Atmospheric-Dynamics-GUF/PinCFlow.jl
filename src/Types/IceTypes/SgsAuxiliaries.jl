struct SgsAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
	saux1::A # ice crystal number concentration n
	saux2::A # ice mixing ratio q
	saux3::A # vapor mixing ratio qv
end

function SgsAuxiliaries(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
)
	(; cloudcover) = namelists.ice

	return SgsAuxiliaries(
		namelists,
		domain,
		subgrid,
		cloudcover,
	)
end

function SgsAuxiliaries(
	namelists::Namelists,
	domain::Domain,	
	subgrid::SubGrid,
	cloudcover::CloudCoverOff,
)
	saux1 = zeros(0, 0, 0)
	saux2 = zeros(0, 0, 0)
	saux3 = zeros(0, 0, 0)

	return SgsAuxiliaries(saux1, saux2, saux3)
end

function SgsAuxiliaries(
	namelists::Namelists,
	domain::Domain,
	subgrid::SubGrid,
	cloudcover::CloudCoverOn,
)

	(; i0, i1, j0, j1, k0, k1) = domain
	(; nscx, nscy, nscz) = namelists.ice
	(; nxnscxx, nynscyy, nznsczz) = subgrid

	saux1 = zeros(nxnscxx, nynscyy, nznsczz)
	saux2 = zeros(nxnscxx, nynscyy, nznsczz)
	saux3 = zeros(nxnscxx, nynscyy, nznsczz)

	return SgsAuxiliaries(saux1, saux2, saux3)

end

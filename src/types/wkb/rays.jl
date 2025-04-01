struct Rays{A <: AbstractArray{<:AbstractFloat, 4}}
  dens::A
  omega::A
  x::A
  y::A
  z::A
  k::A
  l::A
  m::A
  dxray::A
  dyray::A
  dzray::A
  dkray::A
  dlray::A
  dmray::A
  area_xk::A
  area_yl::A
  area_zm::A
end

function Rays(nray_wrk::Integer, nxx::Integer, nyy::Integer, nzz::Integer)
  return Rays([zeros(nray_wrk, nxx, nyy, nzz) for i in 1:17]...)
end

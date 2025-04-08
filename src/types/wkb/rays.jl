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

function get_positions(rijk::CartesianIndex{4}, rays)
  return rays.x[rijk], rays.y[rijk], rays.z[rijk]
end

function get_wavenumbers(rijk::CartesianIndex{4}, rays)
  return rays.k[rijk], rays.l[rijk], rays.m[rijk]
end

function get_physical_extents(rijk::CartesianIndex{4}, rays)
  return rays.dxray[rijk], rays.dyray[rijk], rays.dzray[rijk]
end

function get_spectral_extents(rijk::CartesianIndex{4}, rays)
  return rays.dkray[rijk], rays.dlray[rijk], rays.dmray[rijk]
end

function get_areas(rijk::CartesianIndex{4}, rays)
  return rays.area_xk[rijk], rays.area_yl[rijk], rays.area_zm[rijk]
end

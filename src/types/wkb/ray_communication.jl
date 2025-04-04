struct RayCommunication{
  A <: AbstractVector{<:Symbol},
  B <: AbstractVector{<:AbstractFloat},
}
  fields::A
  send_rays_left::B
  send_rays_right::B
  recv_rays_left::B
  recv_rays_right::B
  send_rays_back::B
  send_rays_forw::B
  recv_rays_back::B
  recv_rays_forw::B
end

function RayCommunication()
  fields = [
    :x,
    :y,
    :z,
    :k,
    :l,
    :m,
    :dxray,
    :dyray,
    :dzray,
    :dkray,
    :dlray,
    :dmray,
    :omega,
    :dens,
  ]
  return RayCommunication(fields, [zeros(length(fields)) for i in 1:8]...)
end

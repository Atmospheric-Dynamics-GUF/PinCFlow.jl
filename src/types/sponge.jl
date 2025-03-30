struct Sponge{A <: AbstractArray{<:AbstractFloat, 3}, B <: AbstractFloat}

  # Damping coefficients.
  alphaunifiedsponge::A
  kr_sp_tfc::A
  kr_sp_w_tfc::A

  # Vertical sponge extent.
  zsponge::B
  dzsponge::B

  # Lateral sponge extent.
  xsponge0::B
  xsponge1::B
  ysponge0::B
  ysponge1::B
  dxsponge::B
  dysponge::B
end

function Sponge(namelists::Namelists, domain::Domain, grid::Grid)

  # Get parameters.
  (; spongeheight) = namelists.sponge
  (; nxx, nyy, nzz) = domain
  (; lx, ly, lz) = grid

  # Initialize the sponge layer coefficients.
  (kr_sp_tfc, kr_sp_w_tfc, alphaunifiedsponge) =
    (zeros(nxx, nyy, nzz) for i in 1:3)

  # Set up the sponge layer.
  dzsponge = spongeheight * (lz[2] - lz[1])
  zsponge = lz[2] - dzsponge
  dxsponge = 0.5 * spongeheight * (lx[2] - lx[1])
  dysponge = 0.5 * spongeheight * (ly[2] - ly[1])
  xsponge0 = lx[1] + dxsponge
  ysponge0 = ly[1] + dysponge
  xsponge1 = lx[2] - dxsponge
  ysponge1 = ly[2] - dysponge

  # Return a Sponge instance.
  return Sponge(
    alphaunifiedsponge,
    kr_sp_tfc,
    kr_sp_w_tfc,
    zsponge,
    dzsponge,
    xsponge0,
    xsponge1,
    ysponge0,
    ysponge1,
    dxsponge,
    dysponge,
  )
end

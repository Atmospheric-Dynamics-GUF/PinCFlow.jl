abstract type AbstractSponge end
struct ExponentialSponge <: AbstractSponge end
struct COSMOSponge <: AbstractSponge end
struct PolynomialSponge <: AbstractSponge end
struct SinusoidalSponge <: AbstractSponge end

struct Sponge{
  A <: OffsetArray{<:AbstractFloat, 3},
  B <: AbstractFloat,
  C <: Integer,
}

  # Damping coefficients.
  alphaunifiedsponge::A
  kr_sp_tfc::A
  kr_sp_w_tfc::A

  # Vertical sponge extent.
  zsponge::B
  dzsponge::B
  ksponge::C

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
  (; spongelayer, spongetype, spongeheight) = namelists.sponge
  (; nz) = domain
  (; lx, ly, lz) = grid

  # Initialize the sponge layer coefficients.
  (kr_sp_tfc, kr_sp_w_tfc, alphaunifiedsponge) = (
    OffsetArray(
      zeros((d.nx + 2, d.ny + 2, d.nz + 2)),
      0:(d.nx + 1),
      0:(d.ny + 1),
      0:(d.nz + 1),
    ) for i in 1:3
  )

  # Set up the sponge layer.
  if unifiedsponge && spongetype == ExponentialSponge()
    ksponge = 1
  else
    ksponge = d.nz - ceil(Integer, spongeheight * d.nz)
  end
  dzsponge = spongeheight * (lz[1] - lz[0])
  zsponge = lz[1] - dzsponge
  dxsponge = 0.5 * spongeheight * (lx[1] - lx[0])
  dysponge = 0.5 * spongeheight * (ly[1] - ly[0])
  xsponge0 = lx[0] + dxsponge
  ysponge0 = ly[0] + dysponge
  xsponge1 = lx[1] - dxsponge
  ysponge1 = ly[1] - dysponge

  # Return a Sponge instance.
  return Sponge(
    alphaunifiedsponge,
    kr_sp_tfc,
    kr_sp_w_tfc,
    xsponge0,
    xsponge1,
    ysponge0,
    ysponge1,
    zsponge,
    ksponge,
    dxsponge,
    dysponge,
    dzsponge,
  )
end

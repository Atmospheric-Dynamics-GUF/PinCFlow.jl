struct MergedRayVolume{A <: AbstractFloat}
  xmin::A
  xmax::A
  ymin::A
  ymax::A
  zmin::A
  zmax::A
  kmin::A
  kmax::A
  lmin::A
  lmax::A
  mmin::A
  mmax::A
  dens::A
end

function MergedRayVolume()
  return MergedRayVolume([0.0 for i in 1:13]...)
end

function MergedRayVolume(
  merge_mode::AbstractMergeMode,
  self::MergedRayVolume,
  x::AbstractFloat,
  dx::AbstractFloat,
  y::AbstractFloat,
  dy::AbstractFloat,
  z::AbstractFloat,
  dz::AbstractFloat,
  k::AbstractFloat,
  dk::AbstractFloat,
  l::AbstractFloat,
  dl::AbstractFloat,
  m::AbstractFloat,
  dm::AbstractFloat,
  axk::AbstractFloat,
  ayl::AbstractFloat,
  azm::AbstractFloat,
  dens::AbstractFloat,
  omir::AbstractFloat,
)
  if self.dens == 0
    return MergedRayVolume(
      x - dx / 2,
      x + dx / 2,
      y - dy / 2,
      y + dy / 2,
      z - dz / 2,
      z + dz / 2,
      k - dk / 2,
      k + dk / 2,
      l - dl / 2,
      l + dl / 2,
      m - dm / 2,
      m + dm / 2,
      merge_wave_action(merge_mode, axk, ayl, azm, dens, omir),
    )
  else
    return MergedRayVolume(
      min(self.xmin, x - dx / 2),
      max(self.xmax, x + dx / 2),
      min(self.ymin, y - dy / 2),
      max(self.ymax, y + dy / 2),
      min(self.zmin, z - dz / 2),
      max(self.zmax, z + dz / 2),
      min(self.kmin, k - dk / 2),
      max(self.kmax, k + dk / 2),
      min(self.lmin, l - dl / 2),
      max(self.kmax, l + dl / 2),
      min(self.mmin, m - dm / 2),
      max(self.mmax, m + dm / 2),
      merge_wave_action(merge_mode, self, axk, ayl, azm, dens, omir),
    )
  end
end
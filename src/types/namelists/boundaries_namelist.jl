abstract type AbstractBoundaries end
struct PeriodicBoundaries <: AbstractBoundaries end
struct SolidWallBoundaries <: AbstractBoundaries end

struct BoundariesNamelist{
  A <: AbstractBoundaries,
  B <: AbstractBoundaries,
  C <: AbstractBoundaries,
}
  xboundaries::A
  yboundaries::B
  zboundaries::C
end

function BoundariesNamelist(;
  xboundaries = PeriodicBoundaries(),
  yboundaries = PeriodicBoundaries(),
  zboundaries = SolidWallBoundaries(),
)
  return BoundariesNamelist(xboundaries, yboundaries, zboundaries)
end

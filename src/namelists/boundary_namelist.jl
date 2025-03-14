abstract type AbstractBoundary end
struct PeriodicBoundary <: AbstractBoundary end
struct SolidWallBoundary <: AbstractBoundary end

struct BoundaryNamelist{
    A<:AbstractBoundary,
    B<:AbstractBoundary,
    C<:AbstractBoundary,
}
    xboundary::A
    yboundary::B
    zboundary::C
end

function BoundaryNamelist(;
  xboundary = PeriodicBoundary(),
  yboundary = PeriodicBoundary(),
  zboundary = SolidWallBoundary(),
)
  return BoundaryNamelist(xboundary, yboundary, zboundary)
end

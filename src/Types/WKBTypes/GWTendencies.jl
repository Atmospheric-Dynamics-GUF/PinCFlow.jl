"""
    GWTendencies

Gravity wave drag and heating tendency fields.

# Fields
- `dudt`: Zonal wind tendency (gravity wave drag)
- `dvdt`: Meridional wind tendency (gravity wave drag)
- `dthetadt`: Potential temperature tendency (gravity wave heating)
"""
struct GWTendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    dudt::A
    dvdt::A
    dthetadt::A
end

function GWTendencies(nxx::Integer, nyy::Integer, nzz::Integer)
    return GWTendencies([zeros(nxx, nyy, nzz) for i in 1:3]...)
end

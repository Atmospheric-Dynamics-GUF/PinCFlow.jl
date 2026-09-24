"""
```julia
Spectrum{A <: AbstractArray{<:AbstractFloat, 4}}
```

Composite type for the spectrum of the wave field.

```julia
Spectrum(
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::Spectrum
```

Construct a `Spectrum` instance, with arrays sized according to the given dimensions.

# Fields

  - `k::A`: Zonal wavenumbers.

  - `l::A`: Meridional wavenumbers.

  - `m::A`: Vertical wavenumbers.

  - `omega::A`: Intrinsic frequencies.

  - `a::A`: Wave-action densities.

# Arguments

  - `wave_modes`: Number of spectral modes per grid cell.

  - `nxx`: Number of subdomain grid points in ``\\hat{x}``-direction.

  - `nyy`: Number of subdomain grid points in ``\\hat{y}``-direction.

  - `nzz`: Number of subdomain grid points in ``\\hat{z}``-direction.
"""
struct Spectrum{A <: AbstractArray{<:AbstractFloat, 4}}
    k::A
    l::A
    m::A
    omega::A
    a::A
end

function Spectrum(
    wave_modes::Integer,
    nxx::Integer,
    nyy::Integer,
    nzz::Integer,
)::Spectrum
    return Spectrum([zeros(wave_modes, nxx, nyy, nzz) for i in 1:5]...)
end

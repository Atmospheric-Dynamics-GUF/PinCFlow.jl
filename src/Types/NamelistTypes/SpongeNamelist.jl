"""
```julia
SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
    E <: NTuple{3, <:AbstractFloat},
}
```

Namelist for parameters describing the sponge.

```julia
SpongeNamelist(;
    spongelayer::Bool = false,
    sponge_uv::Bool = false,
    spongeheight::AbstractFloat = 5.0E-1,
    spongealphaz_dim::AbstractFloat = 1.0E-2,
    spongealphaz_fac::AbstractFloat = 1.0E+0,
    unifiedsponge::Bool = false,
    lateralsponge::Bool = false,
    spongetype::AbstractSponge = PolynomialSponge(),
    spongeorder::Integer = 1,
    cosmosteps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)
```

Construct a `SpongeNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `spongelayer::A`: Switch for enabling Rayleigh-damping in the sponge layer.

  - `sponge_uv::A`: Switch for applying the non-unified sponge to the horizontal wind.

  - `spongeheight::B`: Fractional vertical extent of the sponge.

  - `spongealphaz_dim::B`: Rayleigh-damping coefficient of the unified sponge.

  - `spongealphaz_fac::B`: Rayleigh-damping coefficient of the non-unified sponge, multiplied by the time step.

  - `unifiedsponge::A`: Switch for the unified sponge (provides several profiles and is applied to all prognostic variables).

  - `lateralsponge::A`: Switch for the lateral unified sponge.

  - `spongetype::C`: Profile of the unified sponge.

  - `spongeorder::D`: Order of the polynomial unified sponge.

  - `cosmosteps::D`: Factor by which the time step is mulitplied in the damping coefficient of the COSMO-like unified sponge.

  - `relax_to_mean::A`: Switch for relaxing the wind towards its averages on the terrain-following surfaces. If set to `false`, the wind is relaxed towards `relaxation_wind`.

  - `perturbation_period::B`: Period of an oscillating perturbation on top of `relaxation_wind`.

  - `perturbation_amplitude::B`: Amplitude of an oscillating perturbation on top of `relaxation_wind`.

  - `relaxation_wind::E`: Wind to be obtained through Rayleigh damping in the unified sponge.
"""
struct SpongeNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractSponge,
    D <: Integer,
    E <: NTuple{3, <:AbstractFloat},
}
    spongelayer::A
    sponge_uv::A
    spongeheight::B
    spongealphaz_dim::B
    spongealphaz_fac::B
    unifiedsponge::A
    lateralsponge::A
    spongetype::C
    spongeorder::D
    cosmosteps::D
    relax_to_mean::A
    perturbation_period::B
    perturbation_amplitude::B
    relaxation_wind::E
end

function SpongeNamelist(;
    spongelayer::Bool = false,
    sponge_uv::Bool = false,
    spongeheight::AbstractFloat = 5.0E-1,
    spongealphaz_dim::AbstractFloat = 1.0E-2,
    spongealphaz_fac::AbstractFloat = 1.0E+0,
    unifiedsponge::Bool = false,
    lateralsponge::Bool = false,
    spongetype::AbstractSponge = PolynomialSponge(),
    spongeorder::Integer = 1,
    cosmosteps::Integer = 1,
    relax_to_mean::Bool = true,
    perturbation_period::AbstractFloat = 0.0E+0,
    perturbation_amplitude::AbstractFloat = 0.0E+0,
    relaxation_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
)
    return SpongeNamelist(
        spongelayer,
        sponge_uv,
        spongeheight,
        spongealphaz_dim,
        spongealphaz_fac,
        unifiedsponge,
        lateralsponge,
        spongetype,
        spongeorder,
        cosmosteps,
        relax_to_mean,
        perturbation_period,
        perturbation_amplitude,
        relaxation_wind,
    )
end
